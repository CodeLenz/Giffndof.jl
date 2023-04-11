#
# Heaviside function H(t-t0)
#
function Heaviside(t::T,t0::T) where T
    
    ifelse(t>=t0,one(T),zero(T))

end
