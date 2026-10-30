*! jwselect predict guard 0.7.0 23oct2026
program define jwselect_p
    version 16.1
    if "`e(cmd)'" == "jwselect_sa" {
        di as err "jwselect_sa retains temporary-regressor source results; observation-level predict is not supported"
        di as err "run jwselect first, then use estat cells for treatment-effect summaries"
        exit 198
    }
    di as err "jwselect stores treatment-effect summaries; observation-level predict is not defined"
    di as err "use estat cells, or restore a source model with its own supported prediction method"
    exit 198
end
