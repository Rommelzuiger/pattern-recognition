% EM  Fit a mixture-of-Gaussians using the expectation-maximisation algorithm.
%
%  [L,MU,C,P] = EM (DATA, NMODELS, TYPE, REG, PLOT_MIX, PLOT_RESP, EPS, MAXITER)
%
%  Uses the EM algorithm to fit a mixture of NMODELS models, where each model
%  is of TYPE = 'gauss', 'aligned', or 'circular'.
%  Optional parameters are: 
%    REG       the regularisation constant (dflt. 0)
%    PLOT_MIX  indicating whether progress should be plotted (dflt. 1)
%    PLOT_RESP indicating whether responsibilites should be plotted (dflt. 1)
%    EPS       epsilon, the change in likelihood to stop training (dflt. 1e-5)
%    MAXITER   maximum number of iterations (dlft. 100)
%  Returns the likelihood L, the estimated means MU{} and covariance 
%  matrices C{}, and the prior probabilities P().

function [likelihood, mu, C, mix] = em (data, nmodels, type, regularisation, plot_progress, plot_resp, epsilon, max_iters)

	if (nargin < 8), max_iters = 100;    end;		
	if (nargin < 7), epsilon = 1e-5;     end;		
	if (nargin < 6), plot_resp = 1;      end;
	if (nargin < 5), plot_progress = 1;  end;
	if (nargin < 4), regularisation = 0; end;
	if (nargin < 3), type = 'gauss';     end;
	if (nargin < 2), error ('insufficient number of parameters.'); end;

	gridsize				= 30;

	N = size (data,1);
	A = nmodels;

	psigns = { 'ro' 'ko' 'bo' 'mo' 'co' 'rs' 'ks' 'bs' 'ms' 'cs' 'ro' 'ko' 'bo' 'mo' 'co' 'rs' 'ks' 'bs' 'ms' 'cs' 'ro' 'ko' 'bo' 'mo' 'co' 'rs' 'ks' 'bs' 'ms' 'cs' };
  lsigns = { 'r-' 'k-' 'b-' 'm-' 'c-' 'r-' 'k-' 'b-' 'm-' 'c-' 'r-' 'k-' 'b-' 'm-' 'c-' 'r-' 'k-' 'b-' 'm-' 'c-' 'r-' 'k-' 'b-' 'm-' 'c-' 'r-' 'k-' 'b-' 'm-' 'c-' };
  csigns = { 'r'  'k'  'b'  'm'  'c'  'r'  'k'  'b'  'm'  'c'  'r'  'k'  'b'  'm'  'c'  'r'  'k'  'b'  'm'  'c'  'r'  'k'  'b'  'm'  'c'  'r'  'k'  'b'  'm'  'c'  };
  ssigns = { 'r.' 'k.' 'b.' 'm.' 'c.' 'r.' 'k.' 'b.' 'm.' 'c.' 'r.' 'k-' 'b-' 'm-' 'c-' 'r-' 'k-' 'b-' 'm-' 'c-' 'r-' 'k-' 'b-' 'm-' 'c-' 'r-' 'k-' 'b-' 'm-' 'c-' };

  % -----------------------------------------------------------------------

  t   = +data';
  d   = size(t,1); 		  % Data is 2D...
  q   = 1;              % ...mapped to 1D.

	% For plotting:

	min_data = min(min(data)); max_data = max(max(data)); 

	if (plot_progress == 1)

		if (d ~= 2)
			fprintf (1, 'error: can only plot progress for 2D data\n');
			plot_progress = 0;
		else
			figure(1); clf; figure(2); clf;
			steps = (max_data-min_data)/(gridsize-5);
			min_data = min_data - 2*steps;
			max_data = max_data + 2*steps;
			[xx,yy]  = meshgrid(min_data:steps:max_data);
			tt = [ reshape(xx,1,size(xx,1)*size(xx,2)); ...
	         reshape(yy,1,size(yy,1)*size(yy,2)) ];
		end;

	end;

	tries = 0; retry = 1;

  while ((retry == 1) & (tries <= 5))
  	
    % Initialisation

    mix   = 1/A * ones (1,A);   % 1 dimension  x A models
    R     = rand (N,A);         % N points     x A models
    for i = 1:A
      mu{i}     = rand(d,1)*(max_data-min_data)+min_data;
      C{i}  		= eye(d)*(max_data-min_data);
    end;

    % Initialisation

    done = 0; retry = 0;
  	iter = 0; prev_likelihood = 1.0e20;

    while ((~done) & (~retry))

  		iter = iter + 1;

  		done = 1;

      % Eqn. 21

      pti = zeros(A,N);
      pt  = zeros(1,N);

      for i = 1:A

    		Cdet{i} = det (C{i});

    		if (Cdet{i} < 1e-15)
	  			if (retry == 0) 	
						tries = tries + 1;
	  				fprintf (1, 'ONE OR MORE COVARIANCE MATRICES HAVE BECOME SINGULAR; RETRYING (%d)\n', tries);
					end;
					retry = 1;
    		else
	     		Cinv{i} = inv (C{i});

    			factor = (2*pi)^(-d/2) * (1/sqrt(Cdet{i}));

					if (d == 1)
            pti(i,:)  = factor * exp (-0.5 * ((t - mu{i}*ones(1,N)) .* ...
                                     (Cinv{i} * (t - mu{i}*ones(1,N)))));
          else
           pti(i,:)  = factor * exp (-0.5 * sum ((t - mu{i}*ones(1,N)) .* ...
                                    (Cinv{i} * (t - mu{i}*ones(1,N)))));
					end;
          pt = pt + mix(i) .* pti(i,:);
  			end;
  		end;

  		if (~retry)

    		likelihood = sum(pt);

    		for i = 1:A

    			% Eqn. 21

          R(:,i) 			= ((pti(i,:) .* mix(i)) ./ pt)';
    			sumR(i) 		= sum (R(:,i));

    	    % Eqn. 22

          mix_new(i) 	= sumR(i) / N;
    	
      	  % Eqn. 23

          mu_new{i} 	= sum (((R(:,i) * ones(1,d)) .* t'))' ./ sumR(i);

          C_new{i} 		= ((R(:,i) * ones(1,d))' .* ...
                    		 (t - mu_new{i} * ones(1,N)) * (t - mu_new{i} * ones(1,N))') / ...
    									   (N * mix_new(i));
    		end;

				if (rem (iter,10) == 0)    	
	    		fprintf (1, '[%3d] L: %2.2f (change: %2.2e); sum (P(j|x)) = ', ...
	    							iter, likelihood, abs ((likelihood - prev_likelihood)/likelihood));
	    		for i = 1:A		
	    			fprintf (1, '%2.2f ', sumR(i));
 		   		end;
 		   		fprintf (1, '; P(j) = ');
 		   		for i = 1:A
 		   			fprintf (1, '%2.2f ', mix(i));
 	  	 		end;
 		   		fprintf (1, '\n');
				end;

    		done = (abs ((likelihood - prev_likelihood)/likelihood) < epsilon);

    		prev_likelihood = likelihood;

        % Update

        for i = 1:A
          mix(i)    = mix_new(i);
          mu{i}     = mu_new{i};
          C{i}      = C_new{i} + eye(d) * regularisation;

    			switch (type)
    				case 'gauss'
    				case 'aligned'
    					C{i} = diag(diag(C{i}));
    				case 'circular'
    					C{i} = mean(diag(C{i})) * eye(d);
    				otherwise
    					error ('unknown type: can be gauss, aligned, or circular');
    			end;

        end;

    		% Plot progress

    		if (plot_progress)

					figure (1); clf;
    			scatterd (data); hold on;
    			title (sprintf ('EM; iteration: %d', iter));
    			for i = 1:A
    				plot (mu{i}(1), mu{i}(2), psigns{i});
    				factor = (2*pi)^(-d/2) * (1/sqrt(Cdet{i}));
    	      zz = factor * exp (-0.5 * sum ((tt - mu{i}*ones(1,size(tt,2))) .* ...
                    	         (Cinv{i} * (tt - mu{i}*ones(1,size(tt,2))))));
    				level = 0.25*max(max(zz));
    				contour (xx, yy, reshape(zz,size(xx,1),size(xx,2)), [level level], csigns{i});
    			end;
    			axis ([min_data max_data min_data max_data]); axis square;
    			drawnow;

    		end;

         if (plot_resp)

				  figure (2); clf;
				  axes ('defaultlineLineWidth',2,'defaultlineMarkerSize',16); hold on;
				  for i = 1:A
				    plot (1:N, R(:,i), ssigns{i});
				  end;
				  xlabel ('Point #'); title ('Responsibilities');
				  axis ([0 N 0 1]); set (gca, 'XTickLabel', []);

				end;

  		end;
    end;
	end;

  if (plot_progress)   
   figure(1); clf;
   steps = (max_data-min_data)/(gridsize-5);
   min_data = min_data - 2*steps;
   max_data = max_data + 2*steps;
   [xx,yy]  = meshgrid(min_data:steps:max_data);
   tt = [ reshape(xx,1,size(xx,1)*size(xx,2)); ...
         reshape(yy,1,size(yy,1)*size(yy,2)) ];
   scatterd (data); hold on;
   title (sprintf ('EM; iteration: %d', iter));
   for i = 1:A
      plot (mu{i}(1), mu{i}(2), psigns{i});
      factor = (2*pi)^(-d/2) * (1/sqrt(Cdet{i}));
      zz = factor * exp (-0.5 * sum ((tt - mu{i}*ones(1,size(tt,2))) .* ...
         (Cinv{i} * (tt - mu{i}*ones(1,size(tt,2))))));
      level = 0.25*max(max(zz));
      contour (xx, yy, reshape(zz,size(xx,1),size(xx,2)), [level level], csigns{i});
   end;
   axis ([min_data max_data min_data max_data]); axis square;
   drawnow;
 end;

 if (plot_resp)   
   figure (2); clf;
   axes ('defaultlineLineWidth',2,'defaultlineMarkerSize',16); hold on;
   for i = 1:A
      plot (1:N, R(:,i), ssigns{i});
   end;
   xlabel ('Point #'); title ('Responsibilities');
   axis ([0 N 0 1]); set (gca, 'XTickLabel', []);
 end;

return

