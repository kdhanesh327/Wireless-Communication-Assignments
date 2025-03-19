% Dhanesh Kumar

function main1(typeofmod,M)
tic
     if typeofmod=='MPSK'
         MPSK(str2double(M));
     else typeofmod=='MQAM'
         MQAM(str2double(M));
     end
toc
end