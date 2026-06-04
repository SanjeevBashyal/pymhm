# -*- coding: utf-8 -*-
"""Fortran namelist parsing and template rendering helpers."""

import os
import re


BLOCK_RE = re.compile(r"^\s*&\s*([A-Za-z_][A-Za-z0-9_]*)")
END_RE = re.compile(r"^\s*/\s*(?:!.*)?$")
ASSIGNMENT_RE = re.compile(
    r"^(\s*)([A-Za-z_][A-Za-z0-9_]*(?:\s*\([^)]*\))?)(\s*=\s*)(.*)$")


def canonical_name(name):
    """Return a forgiving key for block and variable matching."""
    key = re.sub(r"[^a-z0-9]", "", str(name or "").lower())
    aliases = {
        "petminus1": "petm1",
        "directrunoff1": "directrunoff1",
    }
    return aliases.get(key, key)


def variable_base(lhs):
    """Return the variable name without any Fortran array index."""
    return str(lhs).split("(", 1)[0].strip()


def split_inline_comment(text):
    """Split a namelist value from an inline comment."""
    quote = None
    for index, char in enumerate(text):
        if char in ("'", '"'):
            if quote == char:
                quote = None
            elif quote is None:
                quote = char
        elif char == "!" and quote is None:
            return text[:index].rstrip(), text[index:]
    return text.rstrip(), ""


def split_csv(text):
    """Split comma-separated namelist values while respecting quotes."""
    parts = []
    current = []
    quote = None
    for char in text:
        if char in ("'", '"'):
            if quote == char:
                quote = None
            elif quote is None:
                quote = char
            current.append(char)
        elif char == "," and quote is None:
            parts.append("".join(current).strip())
            current = []
        else:
            current.append(char)
    final = "".join(current).strip()
    if final or text.endswith(","):
        parts.append(final)
    return parts


def parse_scalar(text):
    """Parse one scalar value from a namelist assignment."""
    value = str(text).strip()
    if len(value) >= 2 and value[0] == value[-1] and value[0] in ("'", '"'):
        return value[1:-1]

    lower = value.lower()
    if lower in (".true.", "true"):
        return True
    if lower in (".false.", "false"):
        return False

    try:
        if re.search(r"[.deE]", value):
            return float(value.replace("D", "E").replace("d", "e"))
        return int(value)
    except ValueError:
        return value


def parse_value(text):
    """Parse a scalar or comma-separated value from a namelist assignment."""
    value_text, _ = split_inline_comment(text)
    parts = split_csv(value_text)
    if len(parts) > 1:
        return [parse_scalar(part) for part in parts]
    if not parts:
        return ""
    return parse_scalar(parts[0])


def template_values(path):
    """Read assignment defaults from an existing namelist template."""
    values = {}
    current_block = None
    if not os.path.exists(path):
        return values

    with open(path, "r", encoding="utf-8") as template_file:
        for line in template_file:
            block_match = BLOCK_RE.match(line)
            if block_match:
                current_block = canonical_name(block_match.group(1))
                values.setdefault(current_block, {})
                continue
            if current_block and END_RE.match(line):
                current_block = None
                continue
            if not current_block:
                continue
            assignment = ASSIGNMENT_RE.match(line)
            if not assignment:
                continue
            name = canonical_name(variable_base(assignment.group(2)))
            values[current_block][name] = parse_value(assignment.group(4))
    return values


def template_block_order(path):
    """Return namelist block names in the order they appear in a template."""
    order = []
    if not os.path.exists(path):
        return order
    with open(path, "r", encoding="utf-8") as template_file:
        for line in template_file:
            block_match = BLOCK_RE.match(line)
            if block_match:
                order.append(canonical_name(block_match.group(1)))
    return order


def format_scalar(value):
    """Format one Python value as a namelist scalar."""
    if isinstance(value, bool):
        return ".true." if value else ".false."
    if isinstance(value, str):
        return '"' + value.replace('"', '\\"') + '"'
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        return f"{value:g}"
    return str(value)


def format_value(value):
    """Format a scalar or list for a namelist assignment."""
    if isinstance(value, (list, tuple)):
        return ", ".join(format_scalar(item) for item in value)
    return format_scalar(value)


def render_template(template_path, values_by_block, include_blocks=None):
    """Render a namelist template by replacing known assignment values."""
    include_keys = None
    if include_blocks is not None:
        include_keys = {canonical_name(name) for name in include_blocks}

    with open(template_path, "r", encoding="utf-8") as template_file:
        lines = template_file.readlines()

    rendered = []
    current_block = None
    current_included = True
    block_buffer = None

    def write_line(line):
        if include_keys is None:
            rendered.append(line)
        elif current_included and block_buffer is not None:
            block_buffer.append(line)

    for line in lines:
        block_match = BLOCK_RE.match(line)
        if block_match:
            current_block = canonical_name(block_match.group(1))
            current_included = (
                include_keys is None or current_block in include_keys)
            if include_keys is None:
                rendered.append(line)
            elif current_included:
                block_buffer = [line]
            else:
                block_buffer = None
            continue

        if current_block and END_RE.match(line):
            write_line(line)
            if include_keys is not None and current_included and block_buffer:
                if rendered and rendered[-1].strip():
                    rendered.append("\n")
                rendered.extend(block_buffer)
                block_buffer = None
            current_block = None
            current_included = True
            continue

        if current_block:
            assignment = ASSIGNMENT_RE.match(line)
            if assignment and current_block in values_by_block:
                key = canonical_name(variable_base(assignment.group(2)))
                if key in values_by_block[current_block]:
                    value_text = format_value(values_by_block[current_block][key])
                    _, comment = split_inline_comment(assignment.group(4))
                    line = (
                        f"{assignment.group(1)}{assignment.group(2)}"
                        f"{assignment.group(3)}{value_text}"
                        f"{(' ' + comment.strip()) if comment else ''}\n")
            write_line(line)
        elif include_keys is None:
            rendered.append(line)

    if not rendered or rendered[-1].endswith("\n"):
        return "".join(rendered)
    return "".join(rendered) + "\n"


def write_rendered_namelist(template_path, output_path, values_by_block,
                            include_blocks=None):
    """Render and save a namelist file."""
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    content = render_template(template_path, values_by_block, include_blocks)
    with open(output_path, "w", encoding="utf-8", newline="\n") as output_file:
        output_file.write(content)
    return output_path
