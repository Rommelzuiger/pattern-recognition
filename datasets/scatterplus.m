% SCATTERPLUS Display 6D scatterplot
% 
%        scatterplus(A)
% 
%  Displays a scatterplot of the first three features of the  dataset
%  A. The next three features are used to scale the ellipses plotted at
%  each point in each of the three directions. The seventh feature, or
%  the class label if present, is used to determine the color.

function scatterplus (a)

	[n,m] = size(a);

	label  = getlab(a)';
	labels = unique(label);

	mn = min(+a); 
	mx = max(+a);
	sc = ones(1,m)./(mx-mn);

	r1_mn = 0.01 / sc(1); r1_range = 0.04 / sc(1);
	r2_mn = 0.01 / sc(2); r2_range = 0.04 / sc(2);
	r3_mn = 0.01 / sc(3); r3_range = 0.04 / sc(3);
  cl_mn = 0;            cl_range =  1;

	if (m <= 2)
		scatter(a);
	elseif (m >= 3)
		for i = labels
			ind  = find(label==i);
			for j = ind
				cl = 0.5; r1 = 1; r2 = 1; r3 = 1;
				if (m >= 4), r1 = r1_mn + r1_range * sc(4)*(+a(j,4)-mn(4)); end;
				if (m >= 5), r2 = r2_mn + r2_range * sc(5)*(+a(j,5)-mn(5)); end;
				if (m >= 6), r3 = r3_mn + r3_range * sc(6)*(+a(j,6)-mn(6)); end;
				if (length(labels) == 1)
					if (m >= 7), cl = cl_mn + cl_range * sc(7)*(+a(j,7)-mn(7)); end;
				else
					cl = cl_mn + cl_range * (i-min(labels))/(max(labels)-min(labels));
				end;
				[x,y,z] = ellipsoid (+a(j,1),+a(j,2),+a(j,3), r1, r2, r3);
				surf(x,y,z,cl*ones(size(z))); hold on;
			end;	
		end;
		shading interp;	
		colormap jet;
		colorbar;
	end;


return

