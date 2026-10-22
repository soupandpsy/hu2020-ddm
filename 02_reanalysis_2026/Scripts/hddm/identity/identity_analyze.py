"""Identity model comparison, diagnostics and predictive checks."""
import os,sys
from pathlib import Path
os.environ["HU2026_TASK"]="identity"
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from common import entry
if __name__=="__main__":entry("identity","analyze")
