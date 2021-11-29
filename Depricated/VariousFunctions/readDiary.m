% This function reads the diary from IPOPT and outputs several results
%
% Author: Antoine Falisse
% Date: 12/19/2018
% 
function [CPU_IPOPT,CPU_NLP,NIter,Cost,Dual_inf,Cons_Viol,Compl,Error_NLP,OptSol] = readDiary(pathDiary)

%% Read Diary and extract information

fid  = fopen(pathDiary,'r');
text = textscan(fid,'%s','Delimiter','');
fid  = fclose(fid);
OptSol = 0;
for i=1:size(text{1},1)
    if strfind(text{1}{i},'Total CPU secs in IPOPT')
        % need to find where the number starts
        eqi = strfind(text{1}{i},'=');
        temp = text{1,1}{i};
        count = 1;
        while temp(eqi+count) == ' '
            count = count + 1;
        end
        eqi1 = eqi+count; 
        CPU_IPOPT=str2double(text{1}{i}(eqi1:end));    
    end
    if strfind(text{1}{i},'Total CPU secs in NLP')
        % need to find where the number starts
        eqi = strfind(text{1}{i},'=');
        temp = text{1,1}{i};
        count = 1;
        while temp(eqi+count) == ' '
            count = count + 1;
        end
        eqi1 = eqi+count;        
        CPU_NLP=str2double(text{1}{i}(eqi1:end));
    end
    if strfind(text{1}{i},'Number of Iterations')
        % need to find where the number starts
        eqi = strfind(text{1}{i},':');
        temp = text{1,1}{i};
        count = 1;
        while temp(eqi+count) == ' '
            count = count + 1;
        end
        eqi1 = eqi+count;        
        NIter=str2double(text{1}{i}(eqi1:end));
    end
    if strfind(text{1}{i},'Objective')
        % need to find where the number starts
        eqi = strfind(text{1}{i},':') + 26; % +26 to get unscaled value
        temp = text{1,1}{i};
        count = 1;
        while temp(eqi+count) == ' '
            count = count + 1;
        end
        eqi1 = eqi+count;        
        Cost=str2double(text{1}{i}(eqi1:end));
    end
    if strfind(text{1}{i},'Dual infeasibility')
        % need to find where the number starts
        eqi = strfind(text{1}{i},':') + 26; % +26 to get unscaled value
        temp = text{1,1}{i};
        count = 1;
        while temp(eqi+count) == ' '
            count = count + 1;
        end
        eqi1 = eqi+count;        
        Dual_inf=str2double(text{1}{i}(eqi1:end));
    end
    if strfind(text{1}{i},'Constraint violation')
        % need to find where the number starts
        eqi = strfind(text{1}{i},':') + 26; % +26 to get unscaled value
        temp = text{1,1}{i};
        count = 1;
        while temp(eqi+count) == ' '
            count = count + 1;
        end
        eqi1 = eqi+count;        
        Cons_Viol=str2double(text{1}{i}(eqi1:end));
    end
    if strfind(text{1}{i},'Complementarity')
        % need to find where the number starts
        eqi = strfind(text{1}{i},':') + 26; % +26 to get unscaled value
        temp = text{1,1}{i};
        count = 1;
        while temp(eqi+count) == ' '
            count = count + 1;
        end
        eqi1 = eqi+count;        
        Compl=str2double(text{1}{i}(eqi1:end));
    end
    if strfind(text{1}{i},'Overall NLP error')
        % need to find where the number starts
        eqi = strfind(text{1}{i},':') + 26; % +26 to get unscaled value
        temp = text{1,1}{i};
        count = 1;
        while temp(eqi+count) == ' '
            count = count + 1;
        end
        eqi1 = eqi+count;        
        Error_NLP=str2double(text{1}{i}(eqi1:end));
    end
    if strfind(text{1}{i},'EXIT: Optimal Solution Found.')
        OptSol=1;
    end
end

end
