function DSMC_surface2(filename,TRI,pX,pY,pZ,hflag)
%converts triangular mesh into DSMC Bird code format file DS3DT.DAT under filename

%open file for reading
fid                     =   fopen(filename,'wt');

% flag error condition if file not opened
if (fid==-1),
    error('Unable to open file for writing.');
else
    %data properties
    [num_TR,c]          =   size(TRI);
    [c,num_PT]          =   size(pX);
    
    %write header
    if hflag
        fprintf(fid,'%.0f\n',1);%no of surfaces %<<<<<,commented out for DS3TM
        fprintf(fid,'%.0f\n',num_TR);%no of points
    end
%     min(TRI)
%     size(pX)
    %write triangle vertices one triangle per line
    for k=1:num_TR
        INDXa      =   TRI(k,1)+1;
        INDXb      =   TRI(k,2)+1;
        INDXc      =   TRI(k,3)+1;
        Xplt       =   [pX(INDXa) pX(INDXb) pX(INDXc)];
        Yplt       =   [pY(INDXa) pY(INDXb) pY(INDXc)];
        Zplt       =   [pZ(INDXa) pZ(INDXb) pZ(INDXc)];
        %this becomes the TRS array in most instances
        fprintf(fid,'%.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f\n',...
                Xplt(1),Yplt(1),Zplt(1),Xplt(2),Yplt(2),Zplt(2),Xplt(3),Yplt(3),Zplt(3));
    end
    
    %close file
    status              =   fclose(fid);

     
    
end

