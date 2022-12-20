# register custom MD units

# system of units
# mass = g/mol = amu
# distance = Angstrom = 10^(-10) m
# time = 10 femtoseconds = 10^(-14) s
# frequency = 10^(14) s^(-1) = 100 THz
# velocity = Angstrom/(10 femtoseconds) = 10^4 m/s
# energy = 100 kJ/mol

# checked: 10/25/22

###########################################################################################

# necessary line for custom units
Unitful.register(@__MODULE__)

# create MD units. @unit symbol unit_label description definition false
@unit m_MD "m_MD" Mass_MD 1u"u" false
@unit d_MD "d_MD" Distance_MD 1u"Å" false
@unit t_MD "t_MD" Time_MD 1u"10fs" false
@unit f_MD "f_MD" Frequency_MD 1u"100THz" false
@unit v_MD "v_MD" Velocity_MD 1u"10^4*m/s" false
@unit e_MD "e_MD" Energy_MD 1u"100kJ/mol" false


###########################################################################################

# create new Unitful type. needed for mapping to labels
@derived_dimension EnergyPerMole Unitful.𝐌*Unitful.𝐋^2/Unitful.𝐓^2/Unitful.𝐍

"""
map Unitful quantity to label. for making graph labels in outputgraph. may change later.
"""
function qtytolabel(qty::Unitful.Quantity)
    types=[Unitful.Mass,Unitful.Length,Unitful.Time,Unitful.Frequency,Unitful.Velocity,EnergyPerMole]
    labels=["Mass","Length","Time","Frequency","Velocity","Energy"]
    idx=findfirst(isa.(qty,types))
    labels[idx]
end
