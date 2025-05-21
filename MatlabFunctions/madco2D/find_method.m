function t=find_method(stinput, argstr, nargs)

curdirMethod=strcat(stinput,'/method');
fidMethod=fopen(curdirMethod,'r');

tline = fgets(fidMethod);                                   % read the first line
while tline ~= -1
    tstr = findstr(tline,argstr);                           % is argstr present?
    if isempty(tstr)                                        % not yet
        tline = fgets(fidMethod);
    else                                                    % argstr is found
        tstr = findstr(tline,cat(2,'##$',argstr,'=( '));    % is it an array?
        if ~isempty(tstr)
            tline = fgets(fidMethod);                       % start reading elements
            tline_final ='';                                % initialize final string
            
            % while next string does not start with ## or $$, i.e. new
            % variable
            while isempty(findstr(tline,'##')) & isempty(findstr(tline,'$$'))
                tline_final = cat(2,tline_final,tline);     % cat final string with a new one
                tline = fgets(fidMethod);                   % read the next line in the file
            end                                             % end while
            t=strread(tline_final,'%f',nargs,'delimiter',' '); % convert string into an array
            fclose(fidMethod);
            return
        else                                                % scalar value
            tstr = findstr(tline,cat(2,'##$',argstr,'='));
            if ~isempty(tstr)
                tstr = cat(2,'##$',argstr,'=');             % string that includes argstr
                strlength = size(tstr,2);                   % get the length its length
                tline(strlength-1:end);                     % read only post value
                t=str2num(tline(strlength+1:end));          % convert it to number
                fclose(fidMethod);
               return
            else
                tline = fgets(fidMethod);
            end
        end
    end
end
t=0;
fclose(fidMethod)