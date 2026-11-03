module ForwardVerification
using ELES, DelimitedFiles
export figure1_configuration

"""The manuscript's C1, single-peaked, ten-layer viscosity profile.

It is an idealized boundary-layer profile, not an observed-data fit. Its
degree mix remains four general quadratics, four affine and two constants.
"""
function figure1_configuration()
    data = readdlm(joinpath(@__DIR__, "..", "data", "figure1_profile.csv"),
                   ',', Float64; skipstart=1)
    edges = [data[:, 2]; data[end, 3]]
    kinds = [(:constant, :linear, :quadratic)[Int(d)+1] for d in data[:, 4]]
    model = mixed_model(edges, kinds)
    values = zeros(model.nparameters)
    for k in eachindex(kinds)
        values[model.left_index[k]] = data[k, 5]
        values[model.right_index[k]] = data[k, 7]
        model.mid_index[k] == 0 || (values[model.mid_index[k]] = data[k, 6])
    end
    @assert model.nparameters == 13
    (model=model, values=values, edges=edges, kinds=kinds, specification=data)
end
end
