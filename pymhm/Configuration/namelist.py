# -*- coding: utf-8 -*-
"""Fortran namelist parsing and template rendering helpers."""
from __future__ import annotations

import os
import re
from typing import Any

NamelistValue = Any
NamelistBlockValues = dict[str, NamelistValue]
NamelistValues = dict[str, NamelistBlockValues]


BLOCK_RE = re.compile(r"^\s*&\s*([A-Za-z_][A-Za-z0-9_]*)")
END_RE = re.compile(r"^\s*/\s*(?:!.*)?$")
ASSIGNMENT_RE = re.compile(
    r"^(\s*)"
    r"([A-Za-z_][A-Za-z0-9_]*(?:\s*\([^)]*\))?"
    r"(?:%[A-Za-z_][A-Za-z0-9_]*)?)"
    r"(\s*=\s*)(.*)$")


def canonical_name(name: object) -> str:
    """Return a forgiving key for block and variable matching."""
    key = re.sub(r"[^a-z0-9]", "", str(name or "").lower())
    aliases = {
        "petminus1": "petm1",
        "directrunoff1": "directrunoff1",
    }
    return aliases.get(key, key)


def variable_base(lhs: object) -> str:
    """Return the variable name without any Fortran array index."""
    return str(lhs).split("(", 1)[0].strip()


def variable_key(lhs: object) -> str:
    """Return a key that preserves numeric Fortran indices when present."""
    full_key = canonical_name(lhs)
    base_key = canonical_name(variable_base(lhs))
    return full_key if full_key != base_key else base_key


def indexed_template_key(lhs: object) -> str | None:
    """Return internal indexed key for a Fortran indexed lhs."""
    base_name = variable_base(lhs)
    base_key = canonical_name(base_name)
    match = re.search(r"\(([^)]*)\)", str(lhs))
    if not match:
        return None

    parts = [part.strip() for part in match.group(1).split(",")]
    numeric_parts = [part for part in parts if part.isdigit()]
    if not numeric_parts:
        return None

    index = numeric_parts[-1]
    suffix = "class" if base_key == "geoparam" else "domain"
    return f"{base_key}__{suffix}{index}"


def split_inline_comment(text: str) -> tuple[str, str]:
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


def split_csv(text: str) -> list[str]:
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


def parse_scalar(text: object) -> NamelistValue:
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


def parse_value(text: str) -> NamelistValue:
    """Parse a scalar or comma-separated value from a namelist assignment."""
    value_text, _ = split_inline_comment(text)
    parts = split_csv(value_text)
    if len(parts) > 1:
        return [parse_scalar(part) for part in parts]
    if not parts:
        return ""
    return parse_scalar(parts[0])


def template_values(path: str) -> NamelistValues:
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
            key = (
                indexed_template_key(assignment.group(2))
                or variable_key(assignment.group(2))
            )
            values[current_block][key] = parse_value(assignment.group(4))
    return values


def template_block_order(path: str) -> list[str]:
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


def format_scalar(value: NamelistValue) -> str:
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


def format_value(value: NamelistValue) -> str:
    """Format a scalar or list for a namelist assignment."""
    if isinstance(value, (list, tuple)):
        flattened = []
        for item in value:
            if isinstance(item, (list, tuple)):
                flattened.extend(item)
            else:
                flattened.append(item)
        return ", ".join(format_scalar(item) for item in flattened)
    return format_scalar(value)


def indexed_assignment_items(
        block_values: NamelistBlockValues,
        base_key: str,
        suffix: str) -> list[tuple[int, NamelistValue]]:
    """Return sorted domain/class assignment values for a base key."""
    prefix = f"{base_key}__{suffix}"
    items = []
    for key, value in block_values.items():
        if not key.startswith(prefix):
            continue
        index_text = key[len(prefix):]
        if not index_text.isdigit():
            continue
        items.append((int(index_text), value))
    return sorted(items)


def geo_param_lhs(template_path: str, index: int) -> str:
    """Return version-specific GeoParam assignment lhs."""
    if "v5.13" in template_path.replace("\\", "/"):
        return f"GeoParam({index},:)"
    return f"GeoParam(:, {index})"


def indexed_lhs(
        template_path: str,
        variable_name: str,
        suffix: str,
        index: int,
        template_lhs: str | None = None) -> str:
    """Return an lhs for a generated indexed assignment."""
    if suffix == "class" and canonical_name(variable_name) == "geoparam":
        return geo_param_lhs(template_path, index)
    if template_lhs:
        match = re.search(r"\(([^)]*)\)", str(template_lhs))
        if match:
            parts = [part.strip() for part in match.group(1).split(",")]
            if len(parts) == 1:
                return f"{variable_name}({index})"
            if (
                    len(parts) == 2
                    and parts[0] == ":"
                    and parts[1].isdigit()):
                return f"{variable_name}(:,{index})"
    return f"{variable_name}(:,{index})"


def render_template(
        template_path: str,
        values_by_block: NamelistValues,
        include_blocks: list[str] | tuple[str, ...] | None = None) -> str:
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

    def write_line(line: str) -> None:
        if include_keys is None:
            rendered.append(line)
        elif current_included and block_buffer is not None:
            block_buffer.append(line)

    emitted_indexed = set()

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
                lhs = assignment.group(2)
                base_name = variable_base(lhs)
                base_key = canonical_name(base_name)
                lhs_key = variable_key(lhs)
                block_values = values_by_block[current_block]
                domain_items = indexed_assignment_items(
                    block_values, base_key, "domain")
                class_items = indexed_assignment_items(
                    block_values, base_key, "class")

                if domain_items or class_items:
                    emit_key = (current_block, base_key)
                    if emit_key not in emitted_indexed:
                        suffix = "domain" if domain_items else "class"
                        indexed_items = domain_items or class_items
                        for index, indexed_value in indexed_items:
                            generated_lhs = indexed_lhs(
                                template_path, base_name, suffix, index, lhs)
                            line_text = (
                                f"{assignment.group(1)}{generated_lhs}"
                                f"{assignment.group(3)}"
                                f"{format_value(indexed_value)}\n")
                            write_line(line_text)
                        emitted_indexed.add(emit_key)
                    continue

                if lhs_key in block_values or base_key in block_values:
                    key = lhs_key if lhs_key in block_values else base_key
                    value_text = format_value(block_values[key])
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


def write_rendered_namelist(
        template_path: str,
        output_path: str,
        values_by_block: NamelistValues,
        include_blocks: list[str] | tuple[str, ...] | None = None) -> str:
    """Render and save a namelist file."""
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    content = render_template(template_path, values_by_block, include_blocks)
    with open(output_path, "w", encoding="utf-8", newline="\n") as output_file:
        output_file.write(content)
    return output_path
