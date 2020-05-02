% Algorithm MSCPSO
%%----------------pseudo_code------------
% Begin
	% Select the size of paprticles for each subswarm
	% Initialize the velocity and position of each particle in the population
	% Evaluate fitness value
	% Find out local best in each sub-swarm and global best optimum in the population
	
	% Do in parallel until the maximum number of iterations has reached
			% {
			% %%
			% Apply fitness adaptive stategy to update the inertia weight
			% Calculate the new velocity of each particle in sub-swarm i,i=1,2,3,4;
			% Update position
			% Evaluate the fitness value in every sub-swarm
			% %%
			% Refresh local best
			% Update global best
			% %%
			% If a guided condition met???
				% Apply diversity stategy
			% endif
			% }
	% Endedo
	
	% Return best solution
%%------------------------------------------
clc;
clear;
close all;
    CostFunction =  @(x) Sphere(x); 



    nVar = 10;
    VarSize = [1 nVar]; 
    VarMin = -10;	% Lower Bound of Decision Variables
    VarMax = 10;    % Upper Bound of Decision Variables
    wmax = 0.9;
    wmin = 0.4;
    w = wmax;
    c1 = 2;
    c2 = 2;
    nPop = 50 ;
    MaxIt = 100;
	alfa1= 1/6;
	alfa2= 1/3;
	alfa3= 1/2;
    ShowIterInfo = true;
    
    MaxVelocity = 0.2*(VarMax-VarMin);
    MinVelocity = -MaxVelocity;
    GlobalBest.Cost = inf;
	 
	% The Particle Template
    empty_particle.Position = [];
    empty_particle.Velocity = [];
    empty_particle.Cost = [];
    empty_particle.Best.Position = [];
    empty_particle.Best.Cost = [];

    % Create Population Array
    slave_particle = repmat(empty_particle, nPop, 4);				%The range of one sub_swarm

   
    emptySlaveGlobalBest.Cost=inf;
    emptySlaveGlobalBest.Position=[];
    SlaveGlobalBest = repmat(emptySlaveGlobalBest,4,1);
    BestCosts = zeros(MaxIt,1);

for j=1:4
	   
        for i=1:nPop
			% Generate Random Solution
			slave_particle(i,j).Position = unifrnd(VarMin, VarMax, VarSize);

			% Initialize Velocity
			slave_particle(i,j).Velocity = zeros(VarSize);

			% Evaluation
			slave_particle(i,j).Cost = CostFunction(slave_particle(i,j).Position);

			% Update the Personal Best
			slave_particle(i,j).Best.Position = slave_particle(i,j).Position;
			slave_particle(i,j).Best.Cost = slave_particle(i,j).Cost;

			% Update Global Best
            if slave_particle(i,j).Best.Cost < SlaveGlobalBest(j).Cost
				SlaveGlobalBest(j) = slave_particle(i,j).Best;
            end
        end
		
        if SlaveGlobalBest(j).Cost < GlobalBest.Cost
            GlobalBest = SlaveGlobalBest(j);
        end
end


 for it=1:MaxIt
    % sub-swarm1 and sub-swarm2
	for j = 1:2
		for i=1:nPop

			% Update Velocity
			slave_particle(i,j).Velocity = w*slave_particle(i,j).Velocity ...
				+ c1*rand(VarSize).*(slave_particle(i,j).Best.Position - slave_particle(i,j).Position) ...
				+ c2*rand(VarSize).*(GlobalBest.Position - slave_particle(i,j).Position);

			% Apply Velocity Limits
			slave_particle(i,j).Velocity = max(slave_particle(i,j).Velocity, MinVelocity);
			slave_particle(i,j).Velocity = min(slave_particle(i,j).Velocity, MaxVelocity);

			% Update Position
			slave_particle(i,j).Position = slave_particle(i,j).Position + slave_particle(i,j).Velocity;

			% Apply Lower and Upper Bound Limits
			slave_particle(i,j).Position = max(slave_particle(i,j).Position, VarMin);
			slave_particle(i,j).Position = min(slave_particle(i,j).Position, VarMax);

			% Evaluation
			slave_particle(i,j).Cost = CostFunction(slave_particle(i,j).Position);

			% Update Personal Best
			if slave_particle(i,j).Cost < slave_particle(i,j).Best.Cost

				slave_particle(i,j).Best.Position = slave_particle(i,j).Position;
				slave_particle(i,j).Best.Cost = slave_particle(i,j).Cost;

				% Update Global Best
				if slave_particle(i,j).Best.Cost < SlaveGlobalBest(j).Cost
					SlaveGlobalBest(j) = slave_particle(i,j).Best;
				end            

			end
		end
		% select global best in sub-swarm1 and sub-swarm2(GlobalBest12)
		% if SlaveGlobalBest(j).Cost < GlobalBest12.Cost
			% GlobalBest12 = SlaveGlobalBest(j);
		% end
	end 

	% sub-swarm3
	for i=1:nPop
			j = 3;
			% Update Velocity
			%slave_particle(i,j).Velocity = w*slave_particle(i,j).Velocity ...
			%	+ c1*rand(VarSize).*(slave_particle(i,j).Best.Position - slave_particle(i,j).Position) ...
			%	+ c2*rand(VarSize).*(SlaveGlobalBest(j).Position - slave_particle(i,j).Position);
            gama1=slave_particle(i,1).Cost;
			gama2=slave_particle(i,2).Cost;
			
			slave_particle(i,j).Velocity = w*(((gama1+gama2)/gama1)*slave_particle(i,1).Velocity + ((gama1+gama2)/gama2)*slave_particle(i,2).Velocity +slave_particle(i,j).Velocity)...
				+ c1*rand(VarSize).*(slave_particle(i,j).Best.Position - slave_particle(i,j).Position) ...
				+ c2*rand(VarSize).*(GlobalBest.Position - slave_particle(i,j).Position);
            
			% Apply Velocity Limits
			slave_particle(i,j).Velocity = max(slave_particle(i,j).Velocity, MinVelocity);
			slave_particle(i,j).Velocity = min(slave_particle(i,j).Velocity, MaxVelocity);

			% Update Position
			slave_particle(i,j).Position =  alfa1*slave_particle(i,j).Position + alfa2*slave_particle(i,j).Best.Position + alfa3*GlobalBest.Position + slave_particle(i,j).Velocity; 
		

			% Apply Lower and Upper Bound Limits
			slave_particle(i,j).Position = max(slave_particle(i,j).Position, VarMin);
			slave_particle(i,j).Position = min(slave_particle(i,j).Position, VarMax);

			% Evaluation
			slave_particle(i,j).Cost = CostFunction(slave_particle(i,j).Position);

			% Update Personal Best
			if slave_particle(i,j).Cost < slave_particle(i,j).Best.Cost

				slave_particle(i,j).Best.Position = slave_particle(i,j).Position;
				slave_particle(i,j).Best.Cost = slave_particle(i,j).Cost;

				% Update Global Best
				if slave_particle(i,j).Best.Cost < SlaveGlobalBest(j).Cost
					SlaveGlobalBest(j) = slave_particle(i,j).Best;
				end            

			end
	end
		

		
	% sub-swarm4
	for i=1:nPop
			j = 4;
			% Update Velocity
			slave_particle(i,j).Velocity = slave_particle(i,1).Velocity + slave_particle(i,2).Velocity - slave_particle(i,3).Velocity;

			% Apply Velocity Limits
			slave_particle(i,j).Velocity = max(slave_particle(i,j).Velocity, MinVelocity);
			slave_particle(i,j).Velocity = min(slave_particle(i,j).Velocity, MaxVelocity);

			% Update Position
			slave_particle(i,j).Position = slave_particle(i,j).Position + slave_particle(i,j).Velocity;

			% Apply Lower and Upper Bound Limits
			slave_particle(i,j).Position = max(slave_particle(i,j).Position, VarMin);
			slave_particle(i,j).Position = min(slave_particle(i,j).Position, VarMax);

			% Evaluation
			slave_particle(i,j).Cost = CostFunction(slave_particle(i,j).Position);

			% Update Personal Best
			if slave_particle(i,j).Cost < slave_particle(i,j).Best.Cost

				slave_particle(i,j).Best.Position = slave_particle(i,j).Position;
				slave_particle(i,j).Best.Cost = slave_particle(i,j).Cost;

				% Update Global Best
				if slave_particle(i,j).Best.Cost < SlaveGlobalBest(j).Cost
					SlaveGlobalBest(j) = slave_particle(i,j).Best;
				end            

			end
	end
		
	% refresh the Global best in the population
	
    for j =1:4
		if(SlaveGlobalBest(j).Cost < GlobalBest.Cost)
			GlobalBest = SlaveGlobalBest(j);
		end
    end
    
    BestCosts(it) = GlobalBest.Cost;
    % Display Iteration Information
    if ShowIterInfo
        disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCosts(it))]);
    end
 end


%% Results

    figure;
    %plot((BestCosts,'LineWidth',2); 
    semilogy(BestCosts,'LineWidth',2);      
    xlabel('Iteration');
    ylabel('Best Cost');
%end %end function        



