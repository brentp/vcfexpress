
local CSQ = {}
CSQ.__index = CSQ

function CSQ.new(fields, header)
    -- header is assumed to be a dictionary with field names as keys and their index as values
    -- e.g. { Allele = 1, Consequence = 2, IMPACT = 3, cDNA_position = 4, AF = 5 }
    local self = setmetatable({}, CSQ)
    
    self.raw_fields = fields
    self.header = header
    self.parsed_fields = nil  -- This will store the parsed fields when needed

    return self
end

function CSQ:parse_fields()
    -- Ensure this function does not trigger __index recursively
    if not rawget(self, "parsed_fields") then
        if type(self.raw_fields) == "string" then
            self.parsed_fields = string.split(self.raw_fields, "|")
        else
            self.parsed_fields = self.raw_fields
        end
    end
end

function CSQ:__tostring()
    -- Ensure fields are parsed before attempting to create a string representation
    if not self.parsed_fields then
        self:parse_fields()
    end

    -- Build the string representation by iterating over the header
    local parts = {}
    for name, index in pairs(self.header) do
        local value = self.parsed_fields[index]
        table.insert(parts, name .. ": " .. (value or "nil"))
    end

    return "{" .. table.concat(parts, ", ") .. "}"
end

function CSQ:__index(key)
    -- Check if the key exists in the metatable (functions like parse_fields)
    local value = rawget(CSQ, key)
    if value then
        return value
    end

    -- Access parsed_fields using rawget to avoid triggering __index
    if not rawget(self, "parsed_fields") then
        self:parse_fields()
    end

    -- Direct lookup for the field index
    local field_index = rawget(self, "header")[key]
    
    if field_index ~= nil then
        local field_value = rawget(self, "parsed_fields")[field_index]
        if NUMBER_FIELDS and NUMBER_FIELDS[key] then
            return tonumber(field_value)
        else
            return field_value
        end
    end

    -- Use rawget to access the value without invoking __index
    return rawget(self, key)
end

CSQS = {}
CSQS.__index = CSQS

function CSQS.new(csq, header)
    if type(csq) == "string" then
        csq = string.split(csq, ",")
    end
    local self = setmetatable({}, CSQS)
    self.csqs = {}
    for _, fields in ipairs(csq) do
        table.insert(self.csqs, CSQ.new(fields, header))
    end
    return self
end

function CSQS:any(f)
    for _, csq in ipairs(self.csqs) do
        if f(csq) then
            return true
        end
    end
    return false
end

function CSQS:__tostring()
    local csqs = {}
    for _, csq in ipairs(self.csqs) do
        table.insert(csqs, tostring(csq))
    end
    return table.concat(csqs, ",\n")
end

function CSQS:all(f)
    for _, csq in ipairs(self.csqs) do
        if not f(csq) then
            return false
        end
    end
    return true
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
        -- results is a table like { 'Consequence', 'Codons', 'Amino_acids', 'Gene', 'SYMBOL', 'Feature', 'EXON', 'PolyPhen', 'SIFT', 'Protein_position', 'BIOTYPE' }
        -- convert here to a dictionary with the index as value
        local header = {}
        for i, v in ipairs(result) do
            header[v] = i
        end
        return header
    else
        return nil, "Format part not found in description."
    end
end

if ({...})[1] == "test" then -- run as luau scripts/csq.lua -a test

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
    for i, v in pairs(parsed_table1) do
        print(i, v)
    end

    print("\nParsed Table 2:")
    for i, v in pairs(parsed_table2) do
        print(i, v)
    end

    print("\nParsed Table 3:")
    for i, v in pairs(parsed_table3) do
        print(i, v)
    end


    desc =
    'Functional annotations: \'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO\'"'

    for i, v in ipairs(parse_description(desc)) do
        print(i, v)
    end



    local header = { Allele = 1, Consequence = 2, IMPACT = 3, cDNA_position = 4, AF = 5 }
    local vals = "A|missense_variant|MODERATE|c.1A>G|1.1e-05"


    local csq = CSQ.new(vals, header)

    assert(csq.Consequence == 'missense_variant', 'Expected: "missense_variant"')
    assert(csq.AF == 1.1e-05)
    assert(csq.IMPACT == 'MODERATE')
    assert(csq.cDNA_position == 'c.1A>G')
    assert(csq.Allele == 'A')
    assert(csq.Consequence == 'missense_variant')

    local vals = "||||"
    local csq = CSQ.new(vals, header)
    assert(csq.Consequence == '')
    assert(csq.IMPACT == '')
    assert(csq.cDNA_position == '')
    assert(csq.AF == nil)
    assert(csq.Allele == '')   

end
