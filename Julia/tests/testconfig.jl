"For debug purposes, mostly"
function default_config_creation()::Tuple{SimPrealloc,SimParams}
    spa = SimParams(
        lambda=632.8e-9,
        Q=4,
        Nq=1024,
        ks=[sind(10.0), sind(20.0), sind(30.0)],
        Lx=10.0e-6,
        Ni=10,
        surf=GaussianSurface(30.0e-9, 100.0e-9),
        rescale=true
    )

    sp = SimPrealloc(spa.Nq, length(spa.ks))

    return sp, spa
end

function default_params_for_surface_testing(surf::T)::Tuple{SimPrealloc,SimParams} where {T<:RandomSurface}
    spa = SimParams(
        lambda=632.8e-9,
        Q=4,
        Nq=1024,
        ks=[sind(10.0), sind(20.0), sind(30.0)],
        Lx=10.0e-6,
        Ni=10,
        surf=surf,
        rescale=true
    )

    sp = SimPrealloc(spa.Nq, length(spa.ks))

    return sp, spa

end