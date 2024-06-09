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
