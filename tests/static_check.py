#!/usr/bin/env python3
"""Packaging/interface checks, NOT a Stata or Mata compiler.

Run from any directory. Failure means a source/package inconsistency; success
does not establish that Stata syntax or the installed upstream interface works.
"""
from pathlib import Path
import ast
import re
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def main():
    subprocess.run([sys.executable, str(ROOT / "build.py"), "--check"], check=True)
    loader = (ROOT / "_jwselect_mata.do").read_text()
    for line in loader.splitlines():
        if line.lstrip().startswith("*"):
            require("/*" not in line, "block-comment opener hidden inside Stata star comment")
    require(loader.count("/*") == loader.count("*/"),
            "unbalanced block comments in generated loader")
    print("PASS: generated loader has balanced block comments and safe Stata header")
    modules = sorted((ROOT / "mata").glob("*.mata"))
    sources = "\n".join(p.read_text() for p in modules)
    # Ignore comments and string literals for conservative structural scans.
    code = re.sub(r'/\*.*?\*/|//[^\n]*|"(?:[^"\\]|\\.)*"', " ", sources, flags=re.S)
    # Mata ++/-- take a variable name, not a subscripted matrix element.
    # Regression guard for the Stata 17 compilation failure reported 6 Oct 2026.
    indexed_increment = r"\][ \t]*(?:\+\+|--)|(?:\+\+|--)[ \t]*[A-Za-z_]\w*[ \t]*\["
    require(not re.search(indexed_increment, code),
            "indexed Mata ++/-- is invalid; use an explicit element assignment")
    print("PASS: no indexed Mata increment/decrement expressions")
    # Project style guard: separate block bodies and closing braces from other
    # statements. The compact `hasempty = 1; break }` failed in Stata/MP 17.
    compact_blocks = [line for line in code.splitlines()
                      if ("{" in line and line.split("{", 1)[1].strip())
                      or ("}" in line and line.strip() != "}")]
    require(not compact_blocks,
            f"Mata block bodies and closing braces need separate lines: {compact_blocks}")
    print("PASS: Mata blocks use separate body and closing-brace lines")
    # Reserved future-feature names have no legitimate syntax role in these modules.
    # Source: Stata [M-2] reswords; actual type/control keywords are excluded here.
    forbidden = set("aggregate array boolean byte case catch const default delegate "
                    "delete double enum explicit export float friend global inline int "
                    "local long namespace new operator polymorphic private protected "
                    "public quad short signed static strL super switch template this "
                    "throw try typedef typename union unsigned using virtual volatile".split())
    collisions = forbidden & set(re.findall(r"\b[A-Za-z_]\w*\b", code))
    require(not collisions, f"reserved Mata names in maintained source: {collisions}")
    print("PASS: no reserved identifier collisions in maintained Mata modules")
    definition = (r"(?m)^(?:void|(?:real|string)\s+(?:scalar|matrix|rowvector|colvector)|"
                  r"struct\s+\w+\s+scalar)\s+(jw_\w+)\s*\(")
    definitions = re.findall(definition, code)
    calls = set(re.findall(r"\b(jw_\w+)\s*\(", code))
    require(len(definitions) == len(set(definitions)), "duplicate maintained Mata function")
    require(calls <= set(definitions), f"undefined package Mata call: {calls-set(definitions)}")
    for opening, closing in (("(", ")"), ("[", "]"), ("{", "}")):
        depth = 0
        for char in code:
            depth += (char == opening) - (char == closing)
            require(depth >= 0, f"unmatched {closing} in maintained Mata source")
        require(depth == 0, f"unmatched {opening} in maintained Mata source")
    print(f"PASS: {len(definitions)} package Mata functions; calls and delimiters consistent")

    main_ado = (ROOT / "jwselect.ado").read_text()
    bridge = (ROOT / "mata" / "jwselect_bridge.mata").read_text()
    writes = set(re.findall(r'jw_put\("(out\w+)"', bridge))
    allocated = set(re.findall(r"\b(out\w+)\b", main_ado))
    require(writes <= allocated, f"unallocated transport names: {writes-allocated}")
    for name in writes:
        require(f"`{name}'" in main_ado, f"bridge output {name} not consumed by ado")
    require('ereturn local estat_cmd jwselect_estat' in main_ado,
            "postestimation dispatch is not registered")
    print(f"PASS: {len(writes)} matrix transport outputs consumed by command")
    require(bridge.rfind("string scalar jw_version()") > bridge.index("void jw_stata()"),
            "completion marker must follow the entire entry point")
    require('capture mata: assert(findexternal("jw_stata()") != NULL)' in main_ado,
            "loader must verify the entry point")
    smoke = (ROOT / "tests" / "smoke.do").read_text()
    require(smoke.index("jwselect, bmat") < smoke.index("do tests/reference_fixture.do"),
            "smoke test must exercise public autoload before direct Mata tests")
    print("PASS: completion marker, entry-point checks, and first-call autoload test")

    for ado in ROOT.glob("*.ado"):
        require(not re.search(r"\b(?:e?return)\s+matrix\s+\w+\s*=\s*[er]\(", ado.read_text()),
                f"return matrix must use an ordinary matrix name in {ado.name}")
    print("PASS: matrix-return statements use named matrices")

    fixture = (ROOT / "tests" / "reference_fixture.do").read_text()
    require("void jw_test_reference()\n{\nstruct jw_model scalar JW_M" in fixture,
            "typed fixture variables must be declared inside a Mata function")
    require("}\njw_test_reference()\nmata drop jw_test_reference()" in fixture,
            "reference fixture must execute and clean up its test function")
    print("PASS: typed reference-fixture variables are scoped inside a test function")

    # Missing() is a scalar count in Mata, not a matrix of missing flags.
    # Both public-interface test drivers must use the common full comparator.
    assertion_helper = (ROOT / "tests" / "release030_assertions.do").read_text()
    require("good = vec(selectindex(vec(A:<.)))" in assertion_helper,
            "matrix comparison must index all finite entries")
    require("vec(A:>=.):==vec(B:>=.)" in assertion_helper,
            "matrix comparison must check missing positions elementwise")
    for filename in ("release030_smoke.do", "release030_concurrency.do"):
        test_source = (ROOT / "tests" / filename).read_text()
        require("do tests/release030_assertions.do" in test_source,
                f"{filename} must use the shared matrix comparator")
        require(not re.search(r"vec\(\s*missing\(|selectindex\(\s*(?:vec\(\s*)?!missing\(",
                              test_source),
                f"{filename} treats Mata missing() as an elementwise mask")
    print("PASS: release and concurrency tests share full missing-aware comparisons")

    partial_fixture = (ROOT / "tests" / "release030_partial_cache.do").read_text()
    require(partial_fixture.count("colshape(") == 3 and "rowshape(" not in partial_fixture,
            "partial-cache fixture needs three-column matrices via colshape()")
    for check in ("rows(C.selectB)==340 & cols(C.selectB)==3",
                  "rows(C.reportB)==340 & cols(C.reportB)==3",
                  "rows(C.EH)==41 & cols(C.EH)==3"):
        require(check in partial_fixture, "partial-cache fixture lacks shape assertion")
    print("PASS: partial-cache fixture constructs and checks the intended matrix shapes")

    manifest = (ROOT / "jwselect.pkg").read_text().splitlines()
    files = [line[2:].strip() for line in manifest if line.startswith(("f ", "F "))]
    require(len(files) == len(set(files)), "duplicate package manifest entry")
    require("F _jwselect_mata.do" in manifest,
            "required .do loader must use uppercase F for ado-directory installation")
    for filename in files:
        require((ROOT / filename).is_file(), f"missing installable file {filename}")
    require("jwselect.ado" in files and "_jwselect_mata.do" in files,
            "package omits command or loader")
    for filename in ("numerical_check.py", "reselect_check.py", "static_check.py"):
        ast.parse((ROOT / "tests" / filename).read_text())
    ast.parse((ROOT / "build.py").read_text())
    print(f"PASS: {len(files)} installable files exist; Python sources parse")
    print("PASS: structural checks only; Stata/Mata compilation and integration NOT executed")


if __name__ == "__main__":
    main()
