from mpmath import mp
import sys

# Max precision
mp.dps = 600 

# INPUT
t_string = '14.134'
t_current = mp.mpf(t_string) 

# Auto-detect starting decimals
if '.' in t_string:
    current_decimals = len(t_string.split('.')[1])
else:
    current_decimals = 0

def get_abs_zeta(t_val):
    return mp.fabs(mp.zeta(mp.mpc('0.5', t_val)))

def strict_truncate(val, decimals):
    # Standardize to fixed point string to chop cleanly
    s = mp.nstr(val, decimals + 5, min_fixed=-mp.inf, max_fixed=mp.inf)
    if '.' in s:
        integer_part, fractional_part = s.split('.')
        truncated_fraction = fractional_part[:decimals]
        return mp.mpf(f"{integer_part}.{truncated_fraction}")
    return val

print(f"{'Lvl':<4} | {'Step (+Bits)':<12} | {'Status':<8} | {'t_address'} | {'Energy'}")
print("-" * 300)

# SAFETY: Global iteration limit
total_steps = 0
max_steps = 200

while total_steps < max_steps:
    total_steps += 1
    
    # 1. CALCULATE STEP SIZE
    step_size = current_decimals
    if step_size > 25: step_size = 25
    if step_size < 1: step_size = 1 
    
    target_decimals = current_decimals + step_size
    
    # 2. PROBE & SLOPE
    probe_depth = target_decimals + 5
    h = mp.power(10, -probe_depth) 
    
    e_now = get_abs_zeta(t_current)
    if e_now == 0:
        print(">>> PERFECT ZERO <<<")
        break
        
    e_plus = mp.fabs(mp.zeta(mp.mpc('0.5', t_current + h)))
    
    diff = e_plus - e_now
    if diff == 0: slope = 1e-300
    else: slope = diff / h
        
    raw_jump = (e_now / slope) 
    
    # 3. THE CLAMP
    movement_cap = mp.power(10, -current_decimals) 
    
    if mp.fabs(raw_jump) > movement_cap:
        # HIT THE WALL
        safe_jump = movement_cap * mp.sign(raw_jump)
        jump_note = "[CAP]"
        next_decimals_state = current_decimals 
        # If we cap, we display at the CURRENT resolution (since we didn't expand)
        display_decimals = current_decimals
    else:
        # FITS IN TRENCH
        safe_jump = raw_jump
        jump_note = "[FIT]"
        next_decimals_state = target_decimals
        # If we fit, we display at the NEW TARGET resolution
        display_decimals = target_decimals

    t_next = t_current - safe_jump
    
    # 4. TRUNCATION
    # We always maintain the 'target' precision in memory to avoid losing progress
    t_current = strict_truncate(t_next, target_decimals)
    
    # OUTPUT
    e_final = get_abs_zeta(t_current)
    
    # FIX: Use 'display_decimals + 5' to ensure we print ALL the new digits
    t_print = mp.nstr(t_current, display_decimals + 5)
    e_print = mp.nstr(e_final, 45)
    
    print(f"{current_decimals:<4} | +{step_size:<12} | {jump_note:<8} | {t_print} | {e_print}")
    
    # Update State
    current_decimals = next_decimals_state

    # 5. EXIT CONDITIONS
    if current_decimals >= mp.dps - 20: 
        print(">>> REACHED MAX PRECISION LIMIT <<<")
        break
        
    if e_final < mp.power(10, -(mp.dps - 20)):
        print(">>> CONVERGED <<<")
        break

print("-" * 300)
final_e = get_abs_zeta(t_current)
print(f"FINAL SINK ENERGY: {mp.nstr(final_e, 150)}")
# Print full precision at the end
print(f"FINAL T VALUE:     {mp.nstr(t_current, current_decimals + 5)}")