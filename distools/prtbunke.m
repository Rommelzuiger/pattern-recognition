% Read Bunke's data to PRTOOLS

function D = prtbunke(fname)

fid = fopen(fname,'r');
if fid == -1,
  error(['File not found: ' fname])
end

ok   = 0;
info = [];
lab  = [];

while ~ok,
  str = fgetl(fid);
  switch str
    case '.CHARACTERBAND TYPE'
      info.datatype = fgets(fid);

    case '.COST FUNCTION'
      str = fgets(fid);
      [info.costfun,info.cost] = strread(str,'%s%f','delimiter',' ');

    case '.CLASS MEMBERSHIP'    % here come labels
      str = fgets(fid);
      lab = strread(str,'%s','delimiter',' ');

    case '.DISTANCE MATRIX'
      D = fscanf(fid,'%f');
      N = sqrt(length(D));
      D = reshape(D,N,N);
      if issym(D,1e-10),
        D = 0.5*(D+D');
      end
      if ~isempty(lab),
        D = dataset(D,lab,'featlab',lab);
        C = classsizes(D);
        % Should it be class frequencies or equal priors?
        if length(C) > 1,
          D = setprior(D,C/sum(C));
        end
        if ~isempty(info),
          D = setuser(D,info);
        end
        p   = findstr(lower(fname),'_norm');
        if ~isempty(p),
          pos = findstr(lower(fname),'/');  % Linux
          if isempty(pos)
            pos = findstr(lower(fname),'\');  % Windows
          end
          if isempty(pos),
            z = 1;
          else
            I = find(pos < p(1));
            z = pos(I(end))+1;
          end
          D = setname(D,['Distance matrix for ' fname(z:p(1)-1)]);
        else
          D = setname(D,'Distance matrix');
        end
      end
      ok = 1;
    otherwise
        ;
  end
  ok = feof(fid);
end
fclose(fid);
return
