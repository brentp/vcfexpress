CSQ = {}
CSQ.__index = CSQ

function CSQ.new(fields, header)
    -- if fields is a string then split on |
    if type(fields) == "string" then
        fields = string.split(fields, "|")
    end
    -- now fill a table with keys from header and values from fields
    local self = setmetatable({}, CSQ)
    for i, field in ipairs(fields) do
        self[header[i]] = field
    end
    return self
end

function CSQ:__tostring()
    local fields = {}
    for k, v in pairs(self) do
        table.insert(fields, k .. "=" .. v)
    end
    return table.concat(fields, ";")
end

NUMBER_FIELDS = { "AF", "AFR_AF", "AMR_AF", "ASN_AF", "EUR_AF", "EAS_AF", "SAS_AF", "AA_AF", "EA_AF",
    "gnomAD_AF", "gnomAD_AFR_AF", "gnomAD_AMR_AF", "gnomAD_ASJ_AF", "gnomAD_EAS_AF", "gnomAD_FIN_AF", "gnomAD_NFE_AF",
    "gnomAD_OTH_AF",
    "gnomAD_SAS_AF", "gnomAD_NFE_NWE_AF", "gnomAD_NFE_SEU_AF", "gnomAD_NFE_SWE_AF", "gnomAD_NFE_ONF_AF",
    "gnomAD_NFE_EST_AF", "gnomAD_NFE_MED_AF", "gnomAD_NFE_SCA_AF",
    "gnomAD_NFE_BAL_AF", "gnomAD_NFE_IBS_AF", "gnomAD_NFE_TSI_AF", "gnomAD_NFE_FOE_AF", "gnomAD_NFE_NWE_AF",
    "gnomAD_NFE_SEU_AF", "gnomAD_NFE_SWE_AF", "gnomAD_NFE_ONF_AF",
    "gnomAD_NFE_EST_AF", "gnomAD_NFE_MED_AF", "gnomAD_NFE_SCA_AF", "gnomAD_NFE_BAL_AF", "gnomAD_NFE_IBS_AF",
    "gnomAD_NFE_TSI_AF", "gnomAD_NFE_FOE_AF", "gnomAD_NFE_NWE_AF",
    "gnomAD_NFE_SEU_AF", "gnomAD_NFE_SWE_AF", "gnomAD_NFE_ONF_AF", "gnomAD_NFE_EST_AF", "gnomAD_NFE_MED_AF",
    "gnomAD_NFE_SCA_AF", "gnomAD_NFE_BAL_AF", "gnomAD_NFE_IBS_AF",
    "gnomAD_NFE_TSI_AF", "gnomAD_NFE_FOE_AF", "gnomAD_NFE_NWE_AF", "gnomAD_NFE_SEU_AF", "gnomAD_NFE_SWE_AF",
    "gnomAD_NFE_ONF_AF", "gnomAD_NFE_EST_AF", "gnomAD_NFE_MED_AF",
    "gnomAD_NFE_SCA_AF", "gnomAD_NFE_BAL_AF", "gnomAD_NFE_IB",
    "MAX_AF", "MAX_AF_POPS", }
-- if the field starts with gnomAD_ also add gnomADe_... and gnomADg_...
add_gnomad = {}
for _, field in ipairs(NUMBER_FIELDS) do
    if string.match(field, "^gnomAD_") then
        print(string.gsub(field, "^gnomAD_", "gnomADe_"))
        add_gnomad[#add_gnomad + 1] = string.gsub(field, "^gnomAD_", "gnomADe_")
        add_gnomad[#add_gnomad + 1] = string.gsub(field, "^gnomAD_", "gnomADg_")
    end
end
-- now add the new fields to the list
for _, field in ipairs(add_gnomad) do
    table.insert(NUMBER_FIELDS, field)
end
-- now convert NUMBER_FIELDS to a set
local number_fields = {}
for _, field in ipairs(NUMBER_FIELDS) do
    number_fields[field] = true
end
NUMBER_FIELDS = number_fields

function CSQ:__index(key)
    if NUMBER_FIELDS[key] then
        return tonumber(self[key])
    else
        return self[key]
    end
end

--[[

header = { "Allele", "Consequence", "IMPACT", "cDNA_position", "AF" }
vals = "A|missense_variant|MODERATE|c.1A>G|1.1e-05"

c = CSQ.new(vals, header)
print(tostring(c))
print(c.Allele)
print(c.AF)

--]]
