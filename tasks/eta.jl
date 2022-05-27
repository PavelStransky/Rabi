include("../Rabi.jl")

const PATH = ""

""" Comparison of the model with parameters mu, nu or eta
    Derivation of mu, nu <--> eta equivalence given in 
    https://onedrive.live.com/view.aspx?resid=7A922556594A97EC%21583&id=documents&wd=target%28Rabi.one%7CB083CE75-23CE-EF48-90B3-4C4AE2145CF6%2FPerturbation%20Jx%7C0D7CFC73-76C5-B348-A2E9-CA2E131C97EA%2F%29
onenote:https://d.docs.live.net/7a922556594a97ec/Documents/Fyzika/Rabi.one#Perturbation%20Jx&section-id={B083CE75-23CE-EF48-90B3-4C4AE2145CF6}&page-id={0D7CFC73-76C5-B348-A2E9-CA2E131C97EA}&object-id={91DB63FE-FCCC-9207-3664-F5C0B544387B}&10
"""

pyplot(size=(1000, 800))

R = 100
λ = 0.75
η = 0.001
μ = 0.001
ν = 0.001

# Initial and final systems
rabi = Rabi(R=R, λ=λ, η=η)
rabi = Rabi(R=R, λ=λ, μ=μ, ν=ν)

op = AllOperators(rabi)
op[end] = :Par => dense(Parity(rabi))
ExpectationValues(rabi, op; saveData=false, maxt=20)

LevelDynamics(Rabi(R=10, δ=0.5, η=2); ylims=(-2, 2))