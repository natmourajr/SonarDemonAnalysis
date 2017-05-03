function [b,a]=le_filtro(arquivo,fs)
%
% Le filtro do MAGS ou do SCC e retorna coeficientes [b,a] para uso no MATLAB
%
% [b,a]=le_filtro(arquivo)
%

global pathname
if nargin < 2
    fs=31250;
end
if nargin<1
    filename=1;
    arquivo = [];
else
    filename=arquivo;
    pname=[];
end
if isempty(arquivo)
    [filename,pname]=uigetfile([pathname '*.dat;*.biq'],'Filtros digitais');
end
if filename
    pathname=pname;
    fid=fopen([pathname filename],'r');
    linha=fgets(fid);
    if ~isempty(findstr(linha,'NB'))
        nb=sscanf(linha(3:end),'%d');
        disp(['Filtro do XSONAR ' [pathname filename]])
        disp([num2str(nb) ' biquads'])
        if nb>0
            linha=fgets(fid);
            for i=1:nb*6
                linha=fgets(fid);
                x(i)=sscanf(linha(4:end),'%f');
            end
            %x=x(2:end);
            sos=reshape(x,6,length(x)/6)';
        else
            sos = [];
        end
    else
        disp('Filtro do MAGS')
        x=fscanf(fid,'%f');%x=x(2:end);
        sos=reshape(x,5,length(x)/5)';
        sos=[sos(:,3:5) ones(size(sos,1),1) -sos(:,1:2)];
        fs=25000;
    end
    [b,a]=sos2tf(sos);
    if isempty(b),b=1;end
    if isempty(a),a=1;end
    if nargout==0
        freqz(b,a,1024,fs);
        subplot(2,1,1)
        title([pathname filename])
        v=axis;axis([v(1:2) -100 0])
        drawnow
    end
    fclose(fid);
end
