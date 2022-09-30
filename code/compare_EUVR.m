function enzTable = compare_EUVR(enzTable)
type = repelem({'*'},height(enzTable),1);
%type = cell(height(enzTable),1);
for i=1:height(enzTable)
    if enzTable.enz_pos(i)>0
        if enzTable.minUsageBio(i)>=0 && enzTable.minUsage(i)>=0
            if enzTable.minUsage(i) == enzTable.minUsageBio(i)
                if enzTable.maxUsage(i) == enzTable.maxUsageBio(i)
                    type(i) = {'Equal'};
                elseif enzTable.maxUsage(i)< enzTable.maxUsageBio(i)
                    type(i) = {'down_subset'};
                else
                    type(i) = {'up_subset'};
                end
            elseif enzTable.minUsage(i) < enzTable.minUsageBio(i)
                if enzTable.maxUsage(i) < enzTable.minUsageBio(i)
                    type(i) = {'down_distinct'};
                else %candidates.maxUsage(i) >= candidates.minUsageBio(i)
                    type(i) = {'down_overlaped'};
                    if enzTable.maxUsage(i) > enzTable.maxUsageBio(i)
                        type(i) = {'overlaped'};
                    end
                end
            else %candidates.minUsage(i) > candidates.minUsageBio(i)
                if enzTable.maxUsageBio(i)  < enzTable.minUsage(i)
                    type(i) = {'up_distinct'};
                else %candidates.maxUsageBio(i)>=candidates.minUsage(i)
                    type(i) = {'up_overlaped'};
                    if enzTable.maxUsageBio(i)>=enzTable.maxUsage(i)
                        type(i) = {'embedded'};
                    end
                end
                
            end
        end
    else
        type(i) = {'none'};
    end
end
enzTable.EUV_comparison = type;
end