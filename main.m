clear 
clc

swarm = 30;
maxiter = 10000;
dim = 30;
runnum=30;


Function_name='F7'; % Name of the test function that can be from F1 to F18


% Load details of the selected benchmark function
[lu,fobj,threshold]=Get_Functions_details(Function_name,dim);


%Algorithms
[bestfit_PSO , gbestfit_PSO ] = PSO(swarm,maxiter,runnum,lu,dim,fobj); 
[bestfit_PSO_L , gbestfit_PSO_L ] = PSO_L(swarm,maxiter,runnum,lu,dim,fobj);
[bestfit_SPSO , gbestfit_SPSO ] = SPSO(swarm,maxiter,runnum,lu,dim,fobj); 
[bestfit_SDPSO , gbestfit_SDPSO ] = SDPSO(swarm,maxiter,runnum,lu,dim,fobj);
[bestfit_AWPSO , gbestfit_AWPSO ] = AWPSO(swarm,maxiter,runnum,lu,dim,fobj);
[bestfit_ASRPPSO , gbestfit_ASRPPSO,sss1,sss2,E_f,div_ASRPPSO ] = ASRPPSO(swarm,maxiter,runnum,lu,dim,fobj); 


%Analysis data
fprintf('PSO of %s: Best= %e ,Mean=%e, Std=%e, Ratio=%f \n', Function_name,min(gbestfit_PSO(maxiter,:)),mean(gbestfit_PSO(maxiter,:)),std(gbestfit_PSO(maxiter,:)),mean(gbestfit_PSO(maxiter,:)<threshold));
fprintf('LIDIW of %s: Best= %e ,Mean=%e, Std=%e, Ratio=%f \n', Function_name,min(gbestfit_PSO_L(maxiter,:)),mean(gbestfit_PSO_L(maxiter,:)),std(gbestfit_PSO_L(maxiter,:)),mean(gbestfit_PSO_L(maxiter,:)<threshold));
fprintf('SPSO of %s: Best= %e ,Mean=%e, Std=%e, Ratio=%f \n', Function_name,min(gbestfit_SPSO(maxiter,:)),mean(gbestfit_SPSO(maxiter,:)),std(gbestfit_SPSO(maxiter,:)),mean(gbestfit_SPSO(maxiter,:)<threshold));
fprintf('SDPSO of %s: Best= %e ,Mean=%e, Std=%e, Ratio=%f \n', Function_name,min(gbestfit_SDPSO(maxiter,:)),mean(gbestfit_SDPSO(maxiter,:)),std(gbestfit_SDPSO(maxiter,:)),mean(gbestfit_SDPSO(maxiter,:)<threshold));
fprintf('AWPSO of %s: Best= %e ,Mean=%e, Std=%e, Ratio=%f \n', Function_name,min(gbestfit_AWPSO(maxiter,:)),mean(gbestfit_AWPSO(maxiter,:)),std(gbestfit_AWPSO(maxiter,:)),mean(gbestfit_AWPSO(maxiter,:)<threshold));
fprintf('ASRPPSO of %s: Best= %e ,Mean=%e, Std=%e, Ratio=%f \n', Function_name,min(gbestfit_ASRPPSO(maxiter,:)),mean(gbestfit_ASRPPSO(maxiter,:)),std(gbestfit_ASRPPSO(maxiter,:)),mean(gbestfit_ASRPPSO(maxiter,:)<threshold));


%Analysis plot
tt=1:maxiter/20:maxiter+1; tt(end)=maxiter;

figure; hold on; axis square; box

plot(tt,log10(bestfit_PSO(tt)'),'r--d','MarkerFaceColor','r');
plot(tt,log10(bestfit_PSO_L(tt)'),'c--s','MarkerFaceColor','c');
plot(tt,log10(bestfit_SPSO(tt)'),'g--+','MarkerFaceColor','g');
plot(tt,log10(bestfit_SDPSO(tt)'),'y-->','MarkerFaceColor','y');
plot(tt,log10(bestfit_AWPSO(tt)'),'m--h','MarkerFaceColor','m');

plot(tt,log10(bestfit_ASRPPSO(tt)'),'k--o','MarkerFaceColor','k');

title(Function_name);
xlabel('Generation number');
ylabel('Mean fitness value');
%legend('PSO','PSO-LIDIW','SPSO','SDPSO','AWPSO','ASRPPSO');
legend('ASRPPSO');
hold off;



