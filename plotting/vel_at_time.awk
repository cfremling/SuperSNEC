# vel_at_time.awk
# Usage:
#   awk -v ttarget=172800 -f vel_at_time.awk Data/Ni_init_frac.dat Data/vel.xg
# Prints: mass, velocity for the block in vel.xg with time closest to ttarget.

FNR==NR {             # First file: Ni_init_frac.dat
    n = FNR           # number of zones
    next
}

# Second file: vel.xg
/^ *\"?Time/ {
    # Finish previous block (if any)
    if (have_block) {
        dt = curr_time - ttarget
        if (dt < 0) dt = -dt
        if (!have_best || dt < best_dt) {
            best_dt   = dt
            have_best = 1
            for (i = 1; i <= k; i++) {
                mass_best[i] = mass_curr[i]
                vel_best[i]  = vel_curr[i]
            }
        }
    }

    # Start new block
    have_block = 1
    k = 0

    # Parse time after '='
    t = 0
    for (i = 1; i <= NF; i++) {
        if ($i == "=") {
            t = $(i+1) + 0
            break
        }
    }
    curr_time = t
    next
}

# Within a block, read up to n zone lines (skip blank lines)
have_block && NF > 0 {
    if (k < n) {
        k++
        mass_curr[k] = $1
        vel_curr[k]  = $2
    }
    next
}

END {
    # Check the last block as well
    if (have_block) {
        dt = curr_time - ttarget
        if (dt < 0) dt = -dt
        if (!have_best || dt < best_dt) {
            best_dt = dt
            for (i = 1; i <= k; i++) {
                mass_best[i] = mass_curr[i]
                vel_best[i]  = vel_curr[i]
            }
        }
    }

    # Output mass, velocity for the best block
    for (i = 1; i <= n; i++) {
        print mass_best[i], vel_best[i]
    }
}

