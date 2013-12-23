%
function a = invsig(a)
a = log(a+realmin) - log(1-a+realmin);
return
