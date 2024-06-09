CSQ = {}
CSQ.__index = CSQ

function CSQ.new(fields, header)
    -- if fields is a string then split on |
    if type(fields) == "string" then
        fields = string.split(fields, "|")
    end
    -- now fill a table with keys from header and values from fields
    local self = setmetatable({}, CSQ)
    for i, h in ipairs(header) do
        self[h] = fields[i]
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
    "MAX_AF", "MAX_AF_POPS", "ALLELE_NUM", "DISTANCE"
}
-- if the field starts with gnomAD_ also add gnomADe_... and gnomADg_...
add_gnomad = {}
for _, field in ipairs(NUMBER_FIELDS) do
    if string.match(field, "^gnomAD_") then
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

function parse_description(description)
    -- Remove quotes and spaces
    description = description:gsub('"', ''):gsub('%s+', '')

    -- Determine the correct split point
    local format_str = nil
    if description:find('Format:') then
        format_str = description:match('Format:(.+)')
    else
        format_str = description:match(':\'(.+)')
    end

    if format_str then
        -- Split by | delimiter and store in a table
        local result = {}
        for value in string.gmatch(format_str, '([^|]+)') do
            -- Remove any trailing characters like ')' or ']'
            value = value:gsub('[%)>\'%]%[]', '')
            -- test if value startswith 'Effect(' and remove it if that's the case.
            if value:find('Effect%(') then
                value = value:gsub('^\'?Effect%(', '')
            end

            table.insert(result, value)
        end
        return result
    else
        return nil, "Format part not found in description."
    end
end

--[[
-- Example usage
local input1 =
'#INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence type as predicted by VEP. Format: Consequence|Codons|Amino_acids|Gene|SYMBOL|Feature|EXON|PolyPhen|SIFT|Protein_position|BIOTYPE">'
local input2 =
'##INFO=<ID=EFF,Number=.,Type=String,Description="Predicted effects for this variant.Format: \'Effect ( Effect_Impact | Functional_Class | Codon_Change | Amino_Acid_change| Amino_Acid_length | Gene_Name | Gene_BioType | Coding | Transcript | Exon [ | ERRORS | WARNINGS ] )\'">'
local input3 =
'##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: \'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO\'">'

local parsed_table1 = parse_description(input1)
local parsed_table2 = parse_description(input2)
local parsed_table3 = parse_description(input3)

print("Parsed Table 1:")
for i, v in ipairs(parsed_table1) do
    print(i, v)
end

print("\nParsed Table 2:")
for i, v in ipairs(parsed_table2) do
    print(i, v)
end

print("\nParsed Table 3:")
for i, v in ipairs(parsed_table3) do
    print(i, v)
end


desc =
'Functional annotations: \'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO\'"'

for i, v in ipairs(parse_description(desc)) do
    print(i, v)
end

--]]
