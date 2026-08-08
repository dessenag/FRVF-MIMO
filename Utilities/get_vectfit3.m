function get_vectfit3()
%% get_vectfit3()
%
% Ensures vectfit3.m (B. Gustavsen, SINTEF Energy Research) is available
% on the MATLAB path before FRVF.m calls it. vectfit3.m is deliberately
% NOT redistributed with this repository: its licence restricts it to
% NON-COMMERCIAL use only, which is incompatible with the GNU GPL v3.0
% that covers the rest of this repository, so every user must obtain it
% directly from the original source and accept its licence terms
% themselves -- this function cannot do that for you.
%
% If vectfit3.m is already on the path, this function returns immediately
% and does nothing. Otherwise it:
%   1. prints the source page, the direct package address, and the
%      licence restriction;
%   2. asks the executor to accept an AUTOMATIC download of that package;
%   3. if accepted, downloads and extracts it, keeping only vectfit3.m;
%   4. if declined, or if the automatic download fails for any reason
%      (network error, moved/changed asset, corrupt archive, ...), falls
%      back to a MANUAL flow: prints step-by-step instructions in the
%      Command Window, THEN opens the source page in the system browser,
%      then asks the executor to locate the downloaded file themselves.
% Either way, the executor is shown the licence terms and the original
% source page before anything is fetched, and must explicitly accept.
%
%% Disclaimer
% This function (get_vectfit3.m) is free software, released by the
% authors of this repository under the GNU General Public License v3.0.
% It grants no rights whatsoever to vectfit3.m itself, which remains
% under its own separate, non-commercial licence.
%
%% Credits
% G. Dessena, Universidad Carlos III de Madrid (gdessena@ing.uc3m.es)
%
%% Changelog
% 2026 - Automatic download (with explicit consent) added as the primary
%        path; interactive manual download kept as a fallback.
%
%% References
% Source and licence for vectfit3.m:
%   https://www.sintef.no/en/software/vector-fitting/downloads/vfit3/
% Citation required by that licence when vectfit3.m is used:
%   [1] B. Gustavsen and A. Semlyen, "Rational approximation of frequency
%       domain responses by Vector Fitting", IEEE Trans. Power Delivery,
%       vol. 14, no. 3, pp. 1052-1061, 1999.
%   [2] B. Gustavsen, "Improving the pole relocating properties of
%       vector fitting", IEEE Trans. Power Delivery, vol. 21, no. 3,
%       pp. 1587-1592, 2006.
%   [3] D. Deschrijver, M. Mrozowski, T. Dhaene, D. De Zutter,
%       "Macromodeling of Multiport Systems Using a Fast Implementation
%       of the Vector Fitting Method", IEEE Microwave and Wireless
%       Components Letters, vol. 18, no. 6, pp. 383-385, 2008.

sourceUrl = 'https://www.sintef.no/en/software/vector-fitting/downloads/vfit3/';
zipUrl    = 'https://www.sintef.no/globalassets/project/vectfit/vfit3.zip';
here = fileparts(mfilename('fullpath'));

if exist('vectfit3', 'file') == 2
    return
end

fprintf('\n===============================================================\n');
fprintf('  vectfit3.m not found on the MATLAB path.\n');
fprintf('  This is third-party code (B. Gustavsen, SINTEF Energy\n');
fprintf('  Research), restricted to NON-COMMERCIAL use only, so it is\n');
fprintf('  not redistributed with this repository.\n');
fprintf('  Original source page (licence terms are stated there):\n');
fprintf('    %s\n', sourceUrl);
fprintf('  Direct package address if you accept below:\n');
fprintf('    %s\n', zipUrl);
fprintf('===============================================================\n\n');

answer = input(['Download and install VFIT3.zip automatically from the address\n' ...
    'above, accepting its non-commercial-use licence? [y/N] '], 's');

installed = false;
if ~isempty(answer) && strncmpi(answer, 'y', 1)
    installed = try_automatic_download(zipUrl, here);
    if ~installed
        fprintf('\nAutomatic download did not succeed. Falling back to manual download.\n\n');
    end
else
    fprintf('Skipping automatic download (not accepted).\n\n');
end

if ~installed
    manual_fallback(sourceUrl, here);
end

rehash path
if exist('vectfit3', 'file') ~= 2
    error('get_vectfit3:notFound', 'vectfit3.m still not found on the path; check %s.', here);
end
fprintf('vectfit3.m installed at %s -- continuing.\n', fullfile(here, 'vectfit3.m'));
end

function ok = try_automatic_download(zipUrl, here)
% Best-effort automatic download + extraction. Returns false (never
% throws) on any failure, so the caller can fall back to the manual flow.
ok = false;
fprintf('Downloading VFIT3.zip ...\n');
tmpZip = [tempname(), '.zip'];

try
    websave(tmpZip, zipUrl, weboptions('Timeout', 20));
catch ME
    fprintf('  Automatic download failed: %s\n', ME.message);
    return
end

try
    extractedList = unzip(tmpZip, here);
    match = extractedList(endsWith(extractedList, 'vectfit3.m'));
    if isempty(match)
        fprintf('  Downloaded archive did not contain vectfit3.m.\n');
    else
        if ~strcmp(match{1}, fullfile(here, 'vectfit3.m'))
            movefile(match{1}, fullfile(here, 'vectfit3.m'));
        end
        % VFIT3.zip also bundles example scripts and papers we don't
        % need in this repository; keep only vectfit3.m.
        extras = setdiff(extractedList, match(1));
        for i = 1:numel(extras)
            if exist(extras{i}, 'file'), delete(extras{i}); end
        end
        ok = true;
        fprintf('  Downloaded and extracted successfully.\n');
    end
catch ME
    fprintf('  Could not extract vectfit3.m from the downloaded archive: %s\n', ME.message);
end

if exist(tmpZip, 'file')
    delete(tmpZip);
end
end

function manual_fallback(sourceUrl, here)
% Interactive manual fallback. Prints full instructions in the Command
% Window BEFORE opening the browser, then asks the executor to locate
% the file they downloaded themselves.
fprintf('=========================== Manual download ============================\n');
fprintf('  1. A browser window will now open at:\n');
fprintf('       %s\n', sourceUrl);
fprintf('  2. Review the licence terms on that page (non-commercial use\n');
fprintf('     only), then click "Download: VFIT3.zip".\n');
fprintf('  3. Come back to this Command Window and answer the prompt below;\n');
fprintf('     you will then be asked to locate the downloaded VFIT3.zip (or\n');
fprintf('     an already-extracted vectfit3.m) via a file picker.\n');
fprintf('==========================================================================\n\n');

try
    web(sourceUrl, '-browser');
catch
    fprintf('(Could not open a browser automatically; please open the link above.)\n');
end

answer = input('Have you downloaded VFIT3.zip and want to locate it now? [y/N] ', 's');
if isempty(answer) || ~strncmpi(answer, 'y', 1)
    error('get_vectfit3:notFound', [ ...
        'vectfit3.m is required by FRVF.m but was not found.\n' ...
        'Download it from %s (accepting its non-commercial licence),\n' ...
        'then place vectfit3.m in:\n  %s'], sourceUrl, here);
end

[pickedFile, pickedPath] = uigetfile({'*.zip;*.m', 'VFIT3.zip or vectfit3.m'}, ...
    'Select the downloaded VFIT3.zip or vectfit3.m');
if isequal(pickedFile, 0)
    error('get_vectfit3:notFound', [ ...
        'No file selected. Download vectfit3.m from %s and place it\n' ...
        'in:\n  %s'], sourceUrl, here);
end

fullPicked = fullfile(pickedPath, pickedFile);
[~, ~, ext] = fileparts(fullPicked);

if strcmpi(ext, '.zip')
    extractedList = unzip(fullPicked, here);
    match = extractedList(endsWith(extractedList, 'vectfit3.m'));
    if isempty(match)
        error('get_vectfit3:notFound', 'vectfit3.m was not found inside %s.', fullPicked);
    end
    if ~strcmp(match{1}, fullfile(here, 'vectfit3.m'))
        movefile(match{1}, fullfile(here, 'vectfit3.m'));
    end
    extras = setdiff(extractedList, match(1));
    for i = 1:numel(extras)
        if exist(extras{i}, 'file'), delete(extras{i}); end
    end
else
    copyfile(fullPicked, fullfile(here, 'vectfit3.m'));
end
end
