# src/deps.jl
module CondaSSL

using CondaPkg

CondaDeps = CondaPkg.read_deps(; file = CondaPkg.cur_deps_file())["deps"]

ENV["JULIA_CONDAPKG_ALLOWED_CHANNELS"] = "conda-forge anaconda"

# 1. Ensure we have the SSL bits
for pkg in ("openssl", "ca-certificates")
    if pkg ∉ keys(CondaDeps)
        CondaPkg.add(pkg)
    end
end

# 2. Point Python at the bundle
prefix   = CondaPkg.envdir()
# most Conda envs put the bundle here:
sslfile1 = joinpath(prefix, "ssl", "cacert.pem")
# fallback (if you also have certifi):
sslfile2 = joinpath(prefix,
    "lib", "python",
    "site‑packages", "certifi", "cacert.pem"
)
ENV["SSL_CERT_FILE"] = isfile(sslfile1) ? sslfile1 : sslfile2

end # module
