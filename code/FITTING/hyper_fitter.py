## This python script generates new AMBER parameters by fitting a reference energy to force field dihedral functions.

## This is a modification on Pablo Arantes' python code from ParametrizANI. (https://github.com/pablo-arantes/ParametrizANI.git) 
## Here are some of the original comments:: 
#;This step fits an empirical energy profile to a reference profile. 
#; We use a Python implementation of the [Rotational Profiler](https://rotprof.lncc.br/index.php) (Rusu *et al.*), which provides an analytical algorithm for computing classical torsional parameters. The fitting uses a **linear least-squares regression** method to determine optimal parameters.
#; Users can select the preferred fitting function to best reproduce the reference conformational energy landscape. Notably, **Rotational Profiler does not constrain the phase shift (δ) to 0° or 180°**, allowing for **asymmetric δ values**. While this improves the fit accuracy, it may reduce parameter transferability across stereoisomers.

#> In the CompBioCat Group we have modified the script and some of the original comments are not completely true. 
#> We have added asymetric phase_shift values, though through a sort of organic form. An option to stick to the 
#> original 0 or 120 phase values will be available in the future.

##! For clarification this code is meant to run only through and with the CompBioCat Group's Protocol. So For obvious reasons GFN2-xTB will not be available. 
#;; This has been created to fit AMBER to AIMNET2 (though more combinations are possible). 
#;; These files will be generated automatically and main error and warnings shall come from the .bash script that calls this python.

import logging
import os
import numpy as np
import sys

# Setup logging
logger = logging.getLogger()
logger.setLevel(logging.INFO)
handler = logging.StreamHandler(sys.stdout)
handler.setFormatter(logging.Formatter('[%(levelname)s] %(message)s'))
logger.addHandler(handler)

# --- 1. Load Data & Manage Angles ---

if len(sys.argv) < 3:
    print("Usage: python our_fitter.py <mm_file> <ref_file> [optional_angle_file]")
    print("Example: python our_fitter.py gaff.dat kcal.dat angles.dat")
    sys.exit(1)

mm_filename = sys.argv[1]
ref_filename = sys.argv[2]

try:
    # Load energy columns (force generic text loading to avoid format errors)
    mm_energies = np.loadtxt(mm_filename)
    ref_energies = np.loadtxt(ref_filename)
    
    # Ensure 1D arrays
    mm_energies = np.atleast_1d(mm_energies).flatten()
    ref_energies = np.atleast_1d(ref_energies).flatten()

except Exception as e:
    logger.error(f"Failed to load energy files: {e}")
    sys.exit(1)

# Check for optional Angle File (3rd Argument)
if len(sys.argv) > 3:
    angle_filename = sys.argv[3]
    try:
        logger.info(f"Reading angles from {angle_filename}...")
        phis = np.loadtxt(angle_filename)
        phis = np.atleast_1d(phis).flatten()
    except Exception as e:
        logger.error(f"Failed to load angle file: {e}")
        sys.exit(1)
else:
    logger.info("No angle file provided. Auto-generating angles (-180 to 180)...")
    N_POINTS = len(mm_energies)
    # Generates N points from -180 to 180 (exclusive of upper bound if endpoint=False)
    # Adjust endpoint=True if your data includes both -180 and +180
    phis = np.linspace(-180, 180, N_POINTS, endpoint=False) 

if len(sys.argv) > 4:
    output_filename = sys.argv[4]
else:
    output_filename = "fitter_results.dat"

# Validation
if len(mm_energies) != len(ref_energies):
    logger.error(f"Mismatch: MM file has {len(mm_energies)} points, Ref file has {len(ref_energies)} points.")
    sys.exit(1)

if len(phis) != len(mm_energies):
    logger.error(f"Mismatch: Angle file has {len(phis)} points, but Energy files have {len(mm_energies)}.")
    sys.exit(1)

# --- 2. Math & Fitting Functions ---

pi = np.pi
SMALL_REAL = 1.0E-8

def fit_dihedral(num_active_coeffs, phis, ref_energies, mm_energies, weights, mult, pshift, phase_offset, is_active):
    nPTS = len(phis)
    N = len(is_active) 
    
    active_indices = [k for k, active in enumerate(is_active) if active]
    
    # Linear Least Squares Matrix Setup: A * c = b
    A = np.zeros((num_active_coeffs, num_active_coeffs))
    b = np.zeros(num_active_coeffs)

    # Pre-calculate basis functions
    basis_vals = np.zeros((N, nPTS))
    for k in range(N):
        if is_active[k]:
            for m in range(nPTS):
                # Convert degrees to radians for calculation
                # Note: The formula corresponds to Ryckaert-Bellemans or similar periodic types
                # Ensure your MD engine supports the 'pshift' (usually just +1 or -1)
                ph = (phis[m] - phase_offset[k]) / 360.0 * 2.0 * pi
                basis_vals[k][m] = (1.0 + pshift[k] * np.cos(mult[k] * ph))

    # Fill A and b
    for i_idx, i_real in enumerate(active_indices):
        for j_idx, j_real in enumerate(active_indices):
            sum_val = np.sum(weights * basis_vals[i_real] * basis_vals[j_real])
            A[i_idx][j_idx] = sum_val
        
        diff = ref_energies - mm_energies
        b[i_idx] = np.sum(weights * basis_vals[i_real] * diff)

    # Solve
    try:
        coeffs = np.linalg.solve(A, b)
    except np.linalg.LinAlgError:
        logger.warning("Singular matrix encountered. Using pseudo-inverse.")
        coeffs = np.dot(np.linalg.pinv(A), b)

    # Calculate fitted curve
    fittedcurve = np.zeros(nPTS)
    for m in range(nPTS):
        correction = 0.0
        for idx, k in enumerate(active_indices):
             correction += coeffs[idx] * basis_vals[k][m]
        fittedcurve[m] = mm_energies[m] + correction

    return fittedcurve, coeffs

def calculate_rmse(predicted, actual):
    return np.sqrt(np.mean((predicted - actual) ** 2))

def calculate_mae(predicted, actual):
    return np.mean(np.abs(predicted - actual))

# --- 3. Main Execution ---

import itertools
from itertools import combinations

# Parameter Definitions (Index Mapping):
# Index: 0  1  2  3  4  5  6
# Mult : 0, 1, 1, 2, 2, 3, 6
mult = [0, 1, 1, 2, 2, 3, 6]
pshift = [1, 1, 1, 1, 1, 1, 1]
#phase guess: > 
phase_offsets = [0, 0, 120, 0, 120, 0, 0]
#possible phases to optimize: >
phase_pool = [0, 45, 90, 120, 180]

def filter_pool(first_phase):
	temp_pool = list(phase_pool)
	if first_phase == 0:
		temp_pool.remove(180)
		temp_pool.remove(0)
	elif first_phase == 180:
		temp_pool.remove(0)
		temp_pool.remove(180)
	else:
		temp_pool.remove(first_phase)
	return temp_pool


values = list(range(7))

run_sequences = [
    list(c)
    for r in range(1, len(values) + 1)
    for c in combinations(values, r)
]

# Baseline adjustment (normalize minima to 0.0)
if len(ref_energies) > 0:
    ref_energies = ref_energies - np.min(ref_energies)
    mm_energies = mm_energies - np.min(mm_energies)

weights = np.ones(len(phis))

# New: List to store results before sorting
results_list = []

logger.setLevel(logging.ERROR)
for run_index, sequence in enumerate(run_sequences, start=1):
    selected_mult = [mult[i] for i in sequence]
    selected_pshift = [pshift[i] for i in sequence]
    # Phase guess
    selected_phase_offset = [phase_offsets[i] for i in sequence]
    selected_is_active = [1] * len(sequence)

    old_error = float('inf')
    convergence = False
    
    while not convergence:
        # Optimization of coefficients
        fitted_curve, fitted_coeff = fit_dihedral(
            len(sequence),
            phis,
            ref_energies,
            mm_energies,
            weights,
            selected_mult,
            selected_pshift,
            selected_phase_offset,
            selected_is_active
        )
    
        # Check convergence
        rmse = calculate_rmse(fitted_curve, ref_energies)
        mae = calculate_mae(fitted_curve, ref_energies)
        error = rmse + mae
        error_diff = abs(error-old_error)
        if error_diff < 1e-6:
             convergence = True
             break
        old_error = error
	
	# Optimize phase shift
        best_phase_set = selected_phase_offset
        best_phase_error = error
        for phases in itertools.product(phase_pool, repeat=len(sequence)):
            valid_set = True
            for i in range(len(sequence)):
                if i != 0:
                    if selected_mult[i] == selected_mult[i-1]:
                        allowed_pool = filter_pool(phases[i-1])
                        if phases[i] not in allowed_pool:
                            valid_set = False
                            break
            if not valid_set:
                continue # next phase_shift angle
            
            hyper_fitted_curve, hyper_fitted_coeff = fit_dihedral(
                len(sequence),
                phis,
                ref_energies,
                mm_energies,
                weights,
                selected_mult,
                selected_pshift,
                list(phases),
                selected_is_active
            )
            hyper_rmse = calculate_rmse(hyper_fitted_curve, ref_energies)
            hyper_mae = calculate_mae(hyper_fitted_curve, ref_energies)
            hyper_error = hyper_rmse + hyper_mae
            
            if hyper_error < best_phase_error:
                best_phase_error = hyper_error
                best_phase_set = list(phases)

        selected_phase_offset = best_phase_set
	

    division_validation = np.abs(hyper_rmse/hyper_mae) - 1
    # Store result
    results_list.append({
        'run_index': run_index,
        'rmse': rmse,
        'rmse/mae': division_validation,
        'sequence': sequence,
	'phases': selected_phase_offset,
        'coefficients': hyper_fitted_coeff
    })

    # Save Output (Files are saved with the original run_index)
    output_suffix = f"run{run_index}"

logger.setLevel(logging.INFO)
   



# Sort the list by error from lowest to highest
sorted_results = sorted(results_list, key=lambda x: x['rmse'])

with open(output_filename, "w") as f:
    f.write("\n--- Results Sorted by RMSE (Lowest to Highest) ---\n\n")
    f.write(f"{'Original Run':<14} {'RMSE (kcal/mol)':<20} {'MAE/RMSE':<20} {'Terms Used (Indices)'} {'Coefficients'}\n")
    f.write("-" * 54 + "\n")
    
    for result in sorted_results:
        f.write(f"{result['run_index']:<14} {result['rmse']:.6f} {result['rmse/mae']:.6f} {result['sequence']} {result['coefficients']}\n")
 
#> print("\nProcessing complete. All fit_coeff_run*.dat and fit_curve_run*.dat files were saved using the original Run Index.")


# --- Find and Print the Smallest Function with RMSE < 0.21 ---

best_result = None

current_best = results_list[0]
for result in results_list:
    if result['rmse'] < 0.51 and result['rmse/mae'] < 0.31:
        best_result = result
        break  # stop at the first (smallest function) that meets the criterion

    if current_best['rmse'] + current_best['rmse/mae'] > result['rmse'] + result['rmse/mae']:
        current_best = result



if best_result:
    terms = [mult[i] for i in best_result['sequence']]
    print("\n--- Smallest Function Below Minimal Error Threshold ---")
    print(f"Original Run: {best_result['run_index']}")
    print(f"RMSE (kcal/mol): {best_result['rmse']:.6f}")
    print(f"RMSE/MAE - 1: {best_result['rmse/mae']:.6f}")
    print(f"Terms Used (Indices): {terms}")
    print(f"Shift Values: {best_result['phases']}")
    print(f"Coefficients (kcal/mol): {best_result['coefficients']}")
else:
    terms = [mult[i] for i in current_best['sequence']]
    print("\nNo function found with RMSE below threshold.")
    print("\n--- Best possible parameters ---")
    print(f"Original Run: {current_best['run_index']}")
    print(f"RMSE (kcal/mol): {current_best['rmse']:.6f}")
    print(f"RMSE/MAE: {current_best['rmse/mae']:.6f}")
    print(f"Terms Used (Indices): {terms}")
    print(f"Shift Values: {current_best['phases']}")
    print(f"Coefficients (kcal/mol): {current_best['coefficients']}")


exit()
