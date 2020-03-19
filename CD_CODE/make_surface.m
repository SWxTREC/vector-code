function [TRI,pX,pY,pZ]=make_surface(filename,plt_flag,ViewDir)
%converts .wrl filename into DSMC .txt format

%open file for reading, its a binary file
fid                     =   fopen(filename,'r');

% flag error condition if file not opened
if (fid==-1)
    error('Unable to open file.');
else
    %scan through header
    ID                  =   0;
    while ID==0 || ID==3
        line            =   fgets(fid);
        ID              =   read_ID(line);
    end
    ID                  =   0;
    
    %scan points
    k                   =   1;
    while ID~=3
        line            =   fgets(fid);
        ID              =   read_ID(line);
        if ID==3
            break
        end
        [PTn,n]         =   read_PTS(line);
        pX(1,k:k+n-1)   =   PTn(1:n,1)';
        pY(1,k:k+n-1)   =   PTn(1:n,2)';
        pZ(1,k:k+n-1)   =   PTn(1:n,3)';
        k               =   k+n;
    end
    
    %scan stuffing
    ID                  =   0;
    while ID==0 || ID==3
        line            =   fgets(fid);
        ID              =   read_ID(line);
    end
    line                =   fgets(fid);
    
    %scan triangles
    ID                  =   0;
    k                   =   1;
    while ID~=3
        line            =   fgets(fid);
        ID              =   read_ID(line);
        if ID==3
            break
        end
        [TRn,n]         =   read_TR(line);
        TRI(k:k+n-1,1)  =   TRn(1:n,1);
        TRI(k:k+n-1,2)  =   TRn(1:n,2);
        TRI(k:k+n-1,3)  =   TRn(1:n,3);
        k               =   k+n;
        
    end
    
    %close file
    status              =   fclose(fid);
    
    if plt_flag
        C          =   ones(size(pZ));
        trisurf(TRI+1,pX,pY,pZ,C)
        axis('equal')
        xlabel('X')
        ylabel('Y')
        zlabel('Z')
        colormap([0.5 0.84 0.9])
        grid off
        view(ViewDir(2)-90,ViewDir(1))
    end
     
end

%------------------internal functions---------------------------------
function ID     =   read_ID(line)
%id's the read line as either
%0-skip
%1-points
%2-triangles
%3-end bracket

ID              =   0;
id_char         =   sscanf(line(1),'%c');

if isequal(id_char,'p') || isequal(id_char,'I')

    id_str      =   sscanf(line(1:5),'%s');
    if isequal(id_str,'point')
        ID      =   1;
    end
    if isequal(id_str,'Index')
        ID      =   2;
    end

end
if isequal(id_char,']')
    ID      =   3;
end
    

function [PT,n] =   read_PTS(line)
%reads in coordinates and returns a N by 3
PT_arr          =   strread(line,'%.7f','delimiter',',');%strread(line,'%.6f','delimiter',',')
[r,c]           =   size(PT_arr);
n               =   floor(r/3);
PT_temp         =   reshape(PT_arr,3,n);
PT              =   PT_temp';



function [TR,n] =   read_TR(line)
%reads in triangle indeces and returns a N by 3
TR_arr          =   strread(line,'%f','delimiter',',');
[r,c]           =   size(TR_arr);
% kk=find(TR_arr==0);
% [szrk,szk]=size(kk);
% if szrk>0
%     szrk
%     szk
%     TR_arr
%     pause
% end
%delete negative ones
indx            =   1;
for k=1:r
    if TR_arr(k,1) ~= -1
        TR(indx,1)  =   TR_arr(k,1);%+1indeces start with one - do this when converting?
        indx        =   indx+1;
    end
end
[r,c]           =   size(TR);
n               =   floor(r/3);
TR_temp         =   reshape(TR,3,n);
TR              =   TR_temp';


