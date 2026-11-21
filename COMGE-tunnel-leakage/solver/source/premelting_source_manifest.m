function manifest = premelting_source_manifest(folder)
%PREMELTING_SOURCE_MANIFEST SHA-256 for archived executed MATLAB sources.
files=dir(fullfile(folder,'*.m')); manifest=struct();
for k=1:numel(files)
    fid=fopen(fullfile(folder,files(k).name),'rb'); assert(fid>=0,'PhaseSwitch:IO','Source cannot be read.');
    cleanup=onCleanup(@()fclose(fid)); bytes=fread(fid,Inf,'*uint8'); clear cleanup
    digest=java.security.MessageDigest.getInstance('SHA-256'); digest.update(bytes);
    value=typecast(digest.digest(),'uint8');
    key=matlab.lang.makeValidName(files(k).name);
    manifest.(key)=lower(reshape(dec2hex(value,2)',1,[]));
end
end
