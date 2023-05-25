# register custom MD units

# system of units
# mass = g/mol = amu
# distance = 1 Å = 10^(-10) m
# time = 10 fs = 10^(-14) s
# frequency = 10^(14) s^(-1) = 100 THz
# velocity = 1 Å/(10 fs) = 10^4 m/s
# energy = 100 kJ/mol

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
@unit e⁻ "e⁻" Charge_MD 1 false # needed for charge graph compatibility
# temp
@unit T_MD "T_MD" Temperature_MD 1u"K" false

###########################################################################################

# create new Unitful type. needed for mapping to labels
@derived_dimension EnergyPerMole Unitful.𝐌*Unitful.𝐋^2/Unitful.𝐓^2/Unitful.𝐍

"""
map Unitful quantity to label. for making graph labels in outputgraph. may change later.
"""
qtytolabel(qty::Unitful.Mass) = "Mass"
qtytolabel(qty::Unitful.Length) = "Length"
qtytolabel(qty::Unitful.Time) = "Time"
qtytolabel(qty::Unitful.Frequency) = "Frequency"
qtytolabel(qty::Unitful.Velocity) = "Velocity"
qtytolabel(qty::EnergyPerMole) = "Energy"
qtytolabel(qty::Unitful.DimensionlessQuantity) = "Charge"
# temp
qtytolabel(qty::Unitful.Temperature) = "Temperature"

###########################################################################################

"""get unit label for use in graphs"""
getgraphunitlabel(qty::Unitful.Mass) = u"u"
getgraphunitlabel(qty::Unitful.Length) = u"Å"
getgraphunitlabel(qty::Unitful.Time) = u"fs"
getgraphunitlabel(qty::Unitful.Frequency) = u"THz"
getgraphunitlabel(qty::Unitful.Velocity) = u"m/s"
getgraphunitlabel(qty::EnergyPerMole) = u"kJ/mol"
getgraphunitlabel(qty::Unitful.DimensionlessQuantity) = u"e⁻"
# temp
getgraphunitlabel(qty::Unitful.Temperature) = u"K"
