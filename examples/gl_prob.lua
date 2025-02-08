-- get the posterior genotype probability.
function genotype_prob(GL, alts)
    -- Return 0.0 immediately if index is -1
    if alts == -1 then
        return 0.0
    end
    if #GL < 3 then
        return 0.0
    end

    -- Compute sum of all likelihoods
    local sum_L = 0
    for i = 1, #GL do
        sum_L = sum_L + 10 ^ GL[i]
    end

    -- Compute probability for the requested index
    return (10 ^ GL[alts + 1]) / sum_L
end

function GT_prob(sample)
    return genotype_prob(sample.GL, sample.alts)
end

--[[
local GL = {-6.3, -2.1, 0}

print("Probability for index 1:", genotype_prob(GL, 1)) --   0.0000057359
print("Probability for index 2:", genotype_prob(GL, 2)) --   0.0909085694
print("Probability for index 3:", genotype_prob(GL, 3)) --   0.9090856945
print("Probability for index -1:", genotype_prob(GL, -1)) -- 0.0
--]]

