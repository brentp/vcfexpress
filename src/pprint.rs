// https://stackoverflow.com/a/42062321
// by Alundaio
// used under: https://creativecommons.org/licenses/by-sa/4.0/
pub const PPRINT: &str = r#"
function pprint(node)
    local cache, stack, output = {},{},{}
    local depth = 1
    local output_str = "{"

    while true do
        local size = 0
        for k,v in pairs(node) do
            size = size + 1
        end

        local cur_index = 1
        for k,v in pairs(node) do
            if (cache[node] == nil) or (cur_index >= cache[node]) then

                if (string.find(output_str,"}",output_str:len())) then
                    output_str = output_str .. ",\n"
                elseif not (string.find(output_str,"\n",output_str:len())) then
                    if output_str:len() > 1 then
                        output_str = output_str .. "\n"
                    end
                end

                -- This is necessary for working with HUGE tables otherwise we run out of memory using concat on huge strings
                table.insert(output,output_str)
                output_str = ""

                local key
                if (type(k) == "number" or type(k) == "boolean") then
                    key = "["..tostring(k).."]"
                else
                    key = "."..tostring(k)
                end

                if (type(v) == "number" or type(v) == "boolean") then
                    output_str = output_str .. string.rep('  ',depth) .. key .. " = "..tostring(v)
                elseif (type(v) == "table") then
                    output_str = output_str .. string.rep('  ',depth) .. key .. " = {\n"
                    table.insert(stack,node)
                    table.insert(stack,v)
                    cache[node] = cur_index+1
                    break
                else
                    output_str = output_str .. string.rep('  ',depth) .. key .. " = '"..tostring(v).."'"
                end

                if (cur_index == size) then
                    -- output_str = output_str .. "\n" .. string.rep('  ',depth-1) .. "}"
                    output_str = output_str .. "}"
                else
                    output_str = output_str .. ","
                end
            else
                -- close the table
                if (cur_index == size) then
                    -- output_str = output_str .. "\n" .. string.rep('  ',depth-1) .. "}"
                    output_str = output_str .. "}"
                end
            end

            cur_index = cur_index + 1
        end

        if (size == 0) then
            -- output_str = output_str .. "\n" .. string.rep('  ',depth-1) .. "}"
            output_str = output_str .. "}"
        end

        if (#stack > 0) then
            node = stack[#stack]
            stack[#stack] = nil
            depth = cache[node] == nil and depth + 1 or depth - 1
        else
            break
        end
    end

    -- This is necessary for working with HUGE tables otherwise we run out of memory using concat on huge strings
    table.insert(output,output_str)
    output_str = table.concat(output)

    print(output_str)
end
        "#;

pub const PRELUDE: &str = r#"
function map(f, t, skip_nil)
    local new_t = {}
    local j = 1
    for i, v in ipairs(t) do
        if v ~= nil or not skip_nil then
            new_t[j] = f(v)
            j = j + 1
        end
    end
    return new_t
end

-- note that this  uses ipairs so only the array portions of the table will be used
function filter(f, t, skip_nil)
    local new_t = {}
    local j = 1
    for i, v in ipairs(t) do
        if v ~= nil or not skip_nil then
            if f(v) then
                new_t[j] = v
                j = j + 1
            end
        end
    end
    return new_t
end

-- note that this  uses ipairs so only the array portions of the table will be used
function all(f, t, skip_nil)
    for i, v in ipairs(t) do
        if (v ~= nil or not skip_nil) and not f(v) then
            return false
        end
    end
    return true
end

-- note that this  uses ipairs so only the array portions of the table will be used
function any(f, t, skip_nil)
    for i, v in ipairs(t) do
        if (v ~= nil or not skip_nil) and f(v) then
            return true
        end
    end
    return false
end

"#;
