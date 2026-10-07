"""Optional, explicitly nonrigorous floating-point matrix-power observations."""
import json
from fractions import Fraction
from pathlib import Path
try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("This optional exploratory script requires NumPy.") from exc

def main():
    path = Path(__file__).resolve().parents[1] / "examples" / "coupled32.json"
    request = json.loads(path.read_text(encoding="utf-8"))
    A = np.array([[float(Fraction(x)) for x in row] for row in request["matrix"]])
    alpha = float(Fraction(request["alpha"]))
    B = np.eye(A.shape[0]) - alpha * A
    print(f"# EXPLORATORY ONLY: NumPy {np.__version__}; rounded SVD norms")
    print("power,spectral_norm_exploratory")
    for m in (1, 8, 9, 10, 11, 13, 30, 31):
        norm = np.linalg.norm(np.linalg.matrix_power(B, m), 2)
        print(f"{m},{norm:.17g}")

if __name__ == "__main__":
    main()
