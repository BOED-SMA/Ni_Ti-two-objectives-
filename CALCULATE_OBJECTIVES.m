function [OBJECTIVES_VALUE, ERROR]=CALCULATE_OBJECTIVES(objectfun,NI_COMPO,VF)
%To apply this method on other project, user should provide their own
%'CALCULATE_OBJECTIVES' function. The function shows here just use SMA
%calculation result stored in the 'SMA_MICROMECHANICAL_MODEL_SIMULATOR.xlsx'

[OUTPUT_VARIABLES, ERROR] = SYSTEMS_MODEL(NI_COMPO,VF);


IMPORTANT_VARIABLES=OUTPUT_VARIABLES(:,13:16);

  if ERROR == 1
        OBJECTIVES_VALUE(1, 1:2) = [666, 666]; % DEFINE 666 IN ORDER TO RECOGNISE WHEN THERE WAS AN ERROR
  else
        OBJECTIVES_VALUE(1, 1:2) = objectfun(IMPORTANT_VARIABLES);
  end



end


function [OUTPUT_VARIABLES, ERROR] = SYSTEMS_MODEL(NI_COMPO,VF)

INPUT = xlsread('SMA_MICROMECHANICAL_MODEL_SIMULATOR.xlsx');  % THIS SIMULATES THE MICROMECHANICAL MODEL
idx = 0;
for m = 1:size(INPUT, 1)
    if INPUT(m, 4) == NI_COMPO && INPUT(m, 9) == VF
        idx = m;
    end
end
if idx == 0
    error('something wrong')
end

OUTPUT_VARIABLES=INPUT(idx, :);  
if INPUT(idx, 13)~=888 && INPUT(idx,13)~=999
ERROR=0;
else
ERROR=1; % any value than zero will do
end

end