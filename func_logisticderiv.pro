function func_logisticderiv,x,par
; logistic curve
; L/(1+exp(-k(x-x0)))
; x0 midpoint
; L curve maximum value
; k logistic growth rate or steepness of the curve
; derivative is f(x)*f(-x)
;l = func_logistic(x,par)
;l2 = func_logistic(-x,par)
;return,l*l2
e = exp(-par[2]*(x-par[1]))
return,par[0]*par[2]*e/(1.0+e)^2
end
