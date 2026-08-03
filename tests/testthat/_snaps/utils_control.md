# cat_control_opt

    Code
      cat_control_opt("foo", "bar")
    Output
      foo: bar

---

    Code
      cat_control_opt("foo", NULL)
    Output
      NULL

# cat_optim_args

    Code
      cat_optim_args(ctrl$global_opt$optimArgs)
    Output
      * optimArgs: 
          * method: "Nelder-Mead"
          * poptim: 0.2
          * pressel: 0.6

