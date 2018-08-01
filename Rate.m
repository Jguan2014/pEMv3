function Rate(real_stateSeq,est_stateSeq)
%compare the accuracy of two state sequence 
pass=0;
for g=1:length(est_stateSeq)

    
    beforestate=real_stateSeq(g); %before pEM
    afterstate=est_stateSeq(g);   %after pEM
    
    if beforestate==afterstate
        pass=pass+1;
    end
end


size=length(est_stateSeq);
disp(['succese rate: ' num2str(pass/size)]);
end