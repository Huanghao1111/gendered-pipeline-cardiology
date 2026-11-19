from pathlib import Path
import json

repo = Path(__file__).resolve().parent
files = sorted((repo / "code").rglob("*.py"))
for path in files:
    source = path.read_text(encoding="utf-8-sig")
    compile(source, str(path), "exec")
result = {
    "status": "PASS_PYTHON_SYNTAX_COMPILATION",
    "method": "Built-in compile() on every analysis Python source file under code/; no bytecode cache written",
    "files_checked": len(files),
    "models_executed": False,
}
print(json.dumps(result, indent=2))
