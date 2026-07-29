function prov = provenance_stamp(cfg, data_file, seed)
%PROVENANCE_STAMP  What produced a result file: code, config, data, seed.
%   Save alongside every result MAT so the run can be identified from the file
%   itself. DATA_FILE and SEED are optional.

    if nargin < 2, data_file = ''; end
    if nargin < 3, seed = []; end

    repo = fileparts(mfilename('fullpath'));

    prov.git_sha   = git_out(repo, 'rev-parse HEAD');
    prov.git_dirty = ~isempty(git_out(repo, 'status --porcelain'));
    prov.git_desc  = git_out(repo, 'log -1 --format=%s');
    prov.cfg_hash  = sha256_bytes(getByteStreamFromArray(sort_fields(cfg)));
    prov.matlab    = version;
    prov.timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    prov.seed      = seed;

    if ~isempty(data_file) && exist(data_file, 'file')
        d = dir(data_file);
        prov.data_file  = data_file;
        prov.data_bytes = d.bytes;
        prov.data_hash  = sha256_file(data_file);
    else
        prov.data_file = data_file; prov.data_bytes = NaN; prov.data_hash = '';
    end
end

function s = sort_fields(s)
%SORT_FIELDS  Order struct fields alphabetically, recursively.
%   getByteStreamFromArray serialises fields in declaration order, so without
%   this two configs that are isequaln hash differently.
    if ~isstruct(s), return; end
    s = orderfields(s);
    f = fieldnames(s);
    for e = 1:numel(s)
        for k = 1:numel(f)
            s(e).(f{k}) = sort_fields(s(e).(f{k}));
        end
    end
end

function s = git_out(repo, args)
    [st, out] = system(sprintf('git --no-pager -C "%s" %s 2>/dev/null', repo, args));
    if st ~= 0, s = ''; else, s = strtrim(out); end
end

function h = sha256_bytes(b)
    md = java.security.MessageDigest.getInstance('SHA-256');
    md.update(uint8(b(:)));
    h = lower(reshape(dec2hex(typecast(md.digest(), 'uint8'), 2).', 1, []));
end

function h = sha256_file(f)
    md = java.security.MessageDigest.getInstance('SHA-256');
    fid = fopen(f, 'r');
    if fid < 0, h = ''; return; end
    while true
        b = fread(fid, 4*1024*1024, '*uint8');
        if isempty(b), break; end
        md.update(b);
    end
    fclose(fid);
    h = lower(reshape(dec2hex(typecast(md.digest(), 'uint8'), 2).', 1, []));
end
