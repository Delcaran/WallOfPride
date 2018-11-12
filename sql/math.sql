USE `ios_db`;

--
-- Function asinh
--
create function asinh (x double)
returns double deterministic
return ln(x+sqrt(1+x*x));

--
-- Function sin2
--
create function sin2 (x double)
returns double deterministic
return ((1-cos(2*x))/2);

--
-- Function atanh
--
create function atanh (x double)
returns double deterministic
return ln( (1+x)/(1-x) )/2;

--
-- Function sinh
--
create function sinh (x double)
returns double deterministic
return (exp(x)-exp(-x))/2;

--
-- Function cosh
--
create function cosh (x double)
returns double deterministic
return (exp(x)+exp(-x))/2;