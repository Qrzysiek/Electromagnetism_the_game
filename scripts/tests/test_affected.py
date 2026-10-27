"""affected.py's view of Rust source (DEC-16): which changes are a function's body (they
select the tests that executed the file) and which are declarations (they select every
test). A misclassified change would let a test that the change affects go unrun."""
import hashlib
import os

import affected


def decl(text):
    return affected.declarations(text)


def code(text):
    return " ".join(affected.tokens(text))


BASE = """
pub const LIMIT: f64 = 2e5;
/// Doc.
pub fn step(x: f64) -> f64 {
    let y = x * 2.0; // a comment
    if y > LIMIT { y } else { x }
}
const fn half(x: u32) -> u32 { x / 2 }
trait Shape { fn area(&self) -> f64; }
impl Shape for f64 { fn area(&self) -> f64 { *self * *self } }
"""


def test_a_body_change_keeps_the_declarations():
    changed = BASE.replace("x * 2.0", "x * 3.0")
    assert code(changed) != code(BASE)
    assert decl(changed) == decl(BASE)


def test_a_method_body_change_keeps_the_declarations():
    changed = BASE.replace("*self * *self", "*self + *self")
    assert decl(changed) == decl(BASE)


def test_a_constant_is_a_declaration():
    assert decl(BASE.replace("2e5", "3e5")) != decl(BASE)


def test_a_const_fn_body_is_a_declaration():
    # It runs at compile time, where coverage cannot see it.
    assert decl(BASE.replace("x / 2 }", "x / 3 }")) != decl(BASE)


def test_a_signature_is_a_declaration():
    assert decl(BASE.replace("pub fn step(x: f64)", "pub fn step(x: f32)")) != decl(BASE)


def test_comments_and_layout_change_nothing():
    changed = BASE.replace("// a comment", "// another comment").replace("/// Doc.", "/// Docs.")
    changed = changed.replace("let y = x * 2.0;", "let  y  =  x*2.0 ;")
    assert code(changed) == code(BASE)
    assert decl(changed) == decl(BASE)


def test_braces_in_strings_and_chars_do_not_end_a_body():
    text = 'fn f() { let s = "}}{"; let r = r#"} "{"#; let c = \'}\'; }\nconst AFTER: u8 = 1;'
    d = decl(text)
    assert "AFTER" in d and "\"}}{\"" not in d
    assert decl(text.replace("const AFTER: u8 = 1;", "const AFTER: u8 = 2;")) != d


def test_lifetimes_are_not_chars():
    text = "fn f<'a>(x: &'a str) -> &'a str { x }\nconst K: u8 = 1;"
    assert affected.tokens(text).count("'a") == 3
    assert "K" in decl(text)


def test_a_function_pointer_type_is_not_a_function():
    text = "type Op = fn(f64) -> f64;\nconst K: u8 = 1;"
    assert decl(text) == code(text)


def test_nested_block_comments_are_skipped():
    text = "/* outer /* inner */ still a comment */ const K: u8 = 1;"
    assert affected.tokens(text) == ["const", "K", ":", "u8", "=", "1", ";"]


def test_the_engines_files_parse():
    # Every engine file tokenizes and has declarations no larger than its code.
    for path, (code_hash, decl_hash) in affected.universe().items():
        text = open(f"{affected.ROOT}/{path}", encoding="utf-8").read()
        assert hashlib.sha256(code(text).encode()).hexdigest() == code_hash, path
        assert len(decl(text)) <= len(code(text)), path


def test_unit_test_modules_are_not_engine_files():
    test_only = affected.test_only_files()
    assert "crates/level/src/tests.rs" in test_only
    assert "crates/level/src/tube/tests.rs" in test_only
    assert not test_only & set(affected.universe())


def test_records_are_shared_by_every_checkout(monkeypatch):
    """The records live in the main checkout, whichever checkout (a worktree too) runs: the
    one whose `.git` is a directory, not a worktree's `.git` file; `EM_AFFECTED_RECORDS`
    names another file."""
    monkeypatch.setenv("EM_AFFECTED_RECORDS", "elsewhere.json")
    assert affected.shared_records() == "elsewhere.json"
    monkeypatch.delenv("EM_AFFECTED_RECORDS")
    path = os.path.normpath(affected.shared_records())
    assert path.endswith(os.path.join("target", "affected", "records.json"))
    main = os.path.dirname(os.path.dirname(os.path.dirname(path)))
    assert os.path.isdir(os.path.join(main, ".git")), main
