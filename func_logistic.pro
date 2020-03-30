function func_logistic,x,par
; logistic curve: par=[L,x0,k]
; L/(1+exp(-k(x-x0)))
; x0 midpoint
; L curve maximum value
; k logistic growth rate or steepness of the curve
return,par[0]/(1.0+exp(-par[2]*(x-par[1])))
end
