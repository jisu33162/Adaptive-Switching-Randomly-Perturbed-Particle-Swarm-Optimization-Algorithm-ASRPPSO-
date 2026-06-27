%ASRPPSO

function [bestfit_ASRPPSO , gbestfit_ASRPPSO,sss1,sss2,E_f,div_ASRPPSO] = ASRPPSO(swarm,maxiter,runnum,lu,dim,fobj)

cputime = 0;


c1 = 2.0;
c2 = 2.0;
w = 0.9;
N_inteval=4;

bestfit = zeros(maxiter,1);
gbestfit = zeros(maxiter, runnum);
div_ASRPPSO = zeros(maxiter,runnum);

for run = 1 : runnum
    
    
    mv = 0.2 * (lu(:, 2) - lu(:, 1));
    Vmin  = repmat(-mv, 1, swarm);
    Vmax  = -Vmin;
    velocity = Vmin + (Vmax - Vmin) .* rand(dim, swarm);
    XRRmin = repmat(lu(:, 1), 1, swarm);
    XRRmax = repmat(lu(:, 2), 1, swarm);
    current_position = XRRmin + (XRRmax - XRRmin) .* rand(dim, swarm);
    
    local_best_position  = current_position;
    current_fitness = fobj(current_position');
    
    local_best_fitness  = current_fitness;
    [global_best_fitness,g] = min(local_best_fitness);
    global_best_position = repmat(local_best_position(:,g), 1, swarm);
    
    R1 = rand(dim, swarm);
    R2 = rand(dim, swarm);
    velocity = w*velocity + c1*(R1.*(local_best_position-current_position)) + c2*(R2.*(global_best_position-current_position));
    current_position = current_position + velocity;
    
    iter = 0 ;
    e1 = zeros(dim, swarm);
    e2 = zeros(dim, swarm);
    g1 = zeros(dim, swarm);
    g2 = zeros(dim, swarm);
    
    tic;
    while iter < maxiter
        iter = iter+1;
        current_fitness = fobj(current_position');
        
        for i = 1 : swarm
            if current_fitness(i) < local_best_fitness(i)
                local_best_fitness(i)  = current_fitness(i);
                local_best_position(:,i) = current_position(:,i)   ;
            end
        end
        local_best_history(:,:,iter) = local_best_position;
        [current_global_best_fitness,g] = min(local_best_fitness);
        if current_global_best_fitness < global_best_fitness
            global_best_fitness = current_global_best_fitness;
            global_best_position = repmat(local_best_position(:,g), 1, swarm);
        end
        global_best_history(:,:,iter) = global_best_position;
        bestfit(iter) = bestfit(iter) + global_best_fitness;
        gbestfit(iter,run)=global_best_fitness;
        
        w = 0.9-(0.9-0.4)*(iter/maxiter);
        
        e1 = local_best_position-current_position;
        e2 = global_best_position-current_position;
        g1 = psotanh(e1,0.035,0.275,0,-1.2);
        g2 = psotanh(e2,0.035,0.275,0,-1.2);
        R1 = rand(dim, swarm);
        R2 = rand(dim, swarm);
        N = 50;
        v1 = 0.3*normrnd(0,0.5,[1 N]);
        v2 = 0.7*normrnd(0,0.5,[1 N]);
        alphatau =randsrc(1,N,[0,1;0.25,0.75]);
        sss1=0;
        sss2=0;
        distance1 = 0;
        distance2 = 0;
        v11=0;
        v22=0;
        d1 = zeros(dim,swarm);
        d2 = zeros(dim,swarm);
        D_matrix=squareform(pdist(current_position'));
        d=(sum(D_matrix))'/swarm;
        E_f=(d(g)-min(d))/(max(d)-min(d));
        
        % Confirm current state
        for i=0:N_inteval-1
            if((i/N_inteval) <= E_f && E_f < ((i+1)/N_inteval))
                ksi=i+1;
            end
        end
        %Randomly Occuring
        for t1 = 1:N
            d1 = alphatau(:,t1)*v1(1,t1);
            v11 = d1 + v11;
            distance1=v11*e1;
            
            sss1=v11;
        end
        
        for t2 = 1:N
            d2=alphatau(:,t2)*v2(1,t2);
            v22 = d2 + v22;
            distance2=v22*e2;
            
            sss2=v22;
        end
        
        switch ksi
            case 1  %convergence
                velocity = w * velocity + g1.*(R1.*(e1))+ ...
                    g2.*(R2.*(e2));
                
            case 2 %exploitation
                velocity = w * velocity + g1.*(R1.*(e1)) + ...
                    g2.*(R2.*(e2))+0.1*distance2+distance1;
                
            case 3 %exploration
                velocity = w * velocity + g1.*(R1.*(e1)) + ...
                    g2.*(R2.*(e2))+...
                    0.1*distance1+distance1;
                
            case 4  %jumping-out
                velocity = w * velocity + g1.*(R1.*(e1)) + ...
                    g2.*(R2.*(e2))+...
                    distance1+distance2 ;
        end
        
        %--------------------------------------------------------------
        for i = 1:swarm
            velocity(:, i) = max(velocity(:,i), Vmin(:,i));
            velocity(:, i) = min(velocity(:,i), Vmax(:,i));
        end
        current_position = current_position + velocity;
        for i = 1:swarm
            current_position(:,i) = max(current_position(:,i), lu(:,1));
            current_position(:,i) = min(current_position(:,i), lu(:,2));
        end
        
                    for j = 1:dim
                xm = current_position(j,:);
                xmp = sum(xm);
                xvp(j,:) = xmp/swarm;
            end
            
            xave = repmat(xvp, 1, dim);
            d2 = 0;
            
            for i = 1:swarm
                d1 = 0;
                for j = 1:dim                    
                    d1 = d1 + norm (current_position(j,i) - xave(j,i), 2);                   
%                     di = norm (current_position(j,:) - xvp(j,:), 2);
                end
                d2 = d2 + d1;
            end
            d3 = d2/swarm;
            
            div_ASRPPSO(iter,run) = d3;
    end
    fprintf('ASRPPSO Run No.%d Done!\n',  run); % Show current results for iteration and number of problem.
    disp(['CPU time: ',num2str(toc)]);
    cputime = cputime + toc;
    
    
end


bestfit_ASRPPSO = bestfit/runnum; % Save best fitness and calculate the average best fitness

gbestfit_ASRPPSO = gbestfit;

end