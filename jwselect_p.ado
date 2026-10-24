*! jwselect predict guard 0.4.0 23oct2026
program define jwselect_p
    version 16.1
    di as err "jwselect stores treatment-effect summaries; observation-level predict is not defined"
    di as err "use estat cells, or restore the original jwdid estimation for prediction"
    exit 198
end
