/* jwselect 0.2.0: shared contracts. All V matrices use sampling covariance. */
struct jw_model {
    real colvector z
    real matrix V
    real matrix postmap
    real matrix hmap
    real colvector groups
    real colvector bases
    real matrix T
    real rowvector w
    real scalar p
    real scalar k
    real scalar nunit
    real scalar gap
}

struct jw_family {
    real matrix masks
    real colvector owners
    real colvector ids
}

struct jw_fit {
    real matrix A
    real matrix targets
    real matrix selected
    real matrix scores
    real matrix cellA
    real matrix cell_selected
    real matrix allA
    real matrix allT
    real matrix allmeta
    /* Actual selection criterion can have a different target from reporting. */
    real matrix selA
    real matrix selT
    real matrix outmap
}
