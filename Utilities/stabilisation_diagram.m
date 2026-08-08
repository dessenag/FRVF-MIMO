function [sysIDReduced, clusterData, modeStab, selectedModes] = ...
    stabilisation_diagram(sysID, baseFileName, plotOption, stabParams, ...
                          psdData, fileHandle, plotTitle, saveFigures, mainFontSize, ...
                          textInterpreter, interactive)
% STABILISATION_DIAGRAM
%   Evaluates the stability of modal parameters (frequency, damping, and
%   mode shapes) across multiple model orders and creates a stabilization
%   diagram. It can also produce cluster diagrams (frequency vs damping).
%
%   SYNTAX:
%   [sysIDReduced, clusterData, modeStab, selectedModes] = stabilisation_diagram(...
%       sysID, baseFileName, plotOption, stabParams, ...
%       psdData, fileHandle, plotTitle, saveFigures, mainFontSize, textInterpreter)
%
%   INPUTS:
%       sysID        - (struct array) ID results for multiple system orders:
%                        .order : model order
%                        .ident : [freq; damp; modeShapes...]
%       baseFileName - (string) Base filename for saving figures
%       plotOption   - (string)  'stab'   => standard Stabilisation Diagrams
%                                 'clust'  => frequency-damping diagram (cluster)
%                                 'off'    => no plots
%       stabParams   - (struct)  Contains stability thresholds:
%                                 .d_max, .d_min
%                                 .epsilon
%                                 .dfr, .dz
%                                 .dMAC
%                                 .NumMAC
%                                Optional fields for plotting:
%                                 .f_min, .f_max
%                                 .cbar_min, .cbar_max
%       psdData      - (struct) PSD data for overlay in plots:
%                        .fr_axis => freq axis (vector)
%                        .data    => PSD amplitude (each column = channel)
%                        .chan    => channel indices to plot
%                        .legend  => optional legend entries (not used here)
%       fileHandle   - (string, optional) Additional filename handle
%                      (e.g., 'run01'). Omit or pass '' to skip.
%       plotTitle    - (string, optional) Custom figure title. If empty,
%                      defaults to 'Stabilisation Diagram'.
%       saveFigures  - (logical, optional) If false, no files are saved.
%                      Defaults to true.
%       mainFontSize - (numeric, optional) Main title font size. Axis titles
%                      will use (mainFontSize - 2), and tick/legend text uses
%                      (mainFontSize - 4). Default = 18 if not provided.
%       textInterpreter - (string, optional) Text interpreter applied to all
%                      text elements in every figure (titles, axis labels,
%                      tick labels, legend, colorbar). Accepted values:
%                        'none'  – plain text (default, safest)
%                        'tex'   – TeX subset (MATLAB default rendering)
%                        'latex' – full LaTeX (requires $ delimiters)
%                      Defaults to the current MATLAB figure default when
%                      omitted or passed as [].
%       interactive  - (logical, optional) If true (default), and a plot
%                      was drawn (plotOption ~= 'off'), opens the "Mode
%                      Selection Panel" and blocks (uiwait) until the
%                      user clicks poles and presses Done -- selectedModes
%                      is then whatever was picked. If false, this GUI
%                      step is skipped entirely and the function returns
%                      immediately with selectedModes empty; the plot (if
%                      any) is still drawn. sysIDReduced -- the stability-
%                      screened poles for every order -- is always fully
%                      computed either way, so callers that don't need
%                      interactive picking can take sysIDReduced(end) as
%                      the stability-screened selection at the highest
%                      tested order without waiting on any GUI.
%
%   OUTPUTS:
%       sysIDReduced  - (struct array) Reduced identification data with only
%                       stable modes for each order.
%       clusterData   - (matrix) [order, freq, damp] for stable poles.
%       modeStab      - (struct array) Stability checks (freq/damping/MAC).
%       selectedModes - (matrix) See full docstring above.
%
%   ------------------------------------------------------------------------
%   Author(s):
%     - Original concept: L. Zanotti Fragonara (Cranfield University)
%     - Extended by:      G. Dessena (Cranfield University)
%     - Revised & Comments: G. Dessena
%       (Assistant Professor, Universidad Carlos III de Madrid)
%     - Interactive GUI selection: G. Dessena
%     - Non-interactive mode (interactive = false, tutorial repository): G. Dessena
%   Licence: GNU General Public License v3.0 (GPL 3.0)
%   Note: when saveFigures = true, the export path in this function uses a
%     Windows-style backslash separator ("ID_RESULTS\..."); on macOS/Linux
%     replace it with fullfile('ID_RESULTS', ...) before enabling exports.
%     saveFigures = false (used throughout tutorial.m) is unaffected.
%   ------------------------------------------------------------------------

%% ============= 1. Handle Optional Inputs =============
if nargin < 6 || isempty(fileHandle)
    fileHandle = '';  % No additional handle
end
if nargin < 7 || isempty(plotTitle)
    plotTitle = 'Stabilisation Diagram';
end
if nargin < 8 || isempty(saveFigures)
    saveFigures = true; % Default: save
end
if nargin < 9 || isempty(mainFontSize)
    mainFontSize = 18;  % Default font size for the main title
end
if nargin < 10 || isempty(textInterpreter)
    textInterpreter = get(0, 'DefaultTextInterpreter');
end
if nargin < 11 || isempty(interactive)
    interactive = true;  % Default: interactive picking (original behaviour)
end
validInterp = {'none', 'tex', 'latex'};
if ~any(strcmpi(textInterpreter, validInterp))
    warning('stabilisation_diagram:badInterpreter', ...
        'textInterpreter ''%s'' not recognised; using ''none''.', textInterpreter);
    textInterpreter = 'none';
end

% Derived font sizes
axisFontSize   = mainFontSize - 2; % For x/y axis labels
tickFontSize   = mainFontSize - 4; % For tick labels
legendFontSize = mainFontSize - 6; % For legend

% Ensure plotting limits exist
if ~isfield(stabParams, 'f_min'),     stabParams.f_min    = 0;    end
if ~isfield(stabParams, 'f_max'),     stabParams.f_max    = 500;  end
if ~isfield(stabParams, 'cbar_min'),  stabParams.cbar_min = 0.01; end
if ~isfield(stabParams, 'cbar_max'),  stabParams.cbar_max = 0.10; end

% textInterpreter is now the single source of truth for all text rendering
% in this function's figures (set via nargin block above).

% Initialise outputs / shared GUI state unconditionally so that
% MATLAB's nested-function static scoping creates the shared binding
% for ALL variables used in callbacks, regardless of plotOption.
selectedModes = struct('id', [], 'order', []);  % struct output
selMat        = zeros(0, 3);  % internal GUI: [order, freq_Hz, damp]

%% ============= 2. Remove Poles with Invalid Damping =============
for iOrder = 1:numel(sysID)
    tooLarge = sysID(iOrder).ident(2,:) > stabParams.d_max;
    sysID(iOrder).ident(:, tooLarge) = [];

    tooSmall = sysID(iOrder).ident(2,:) < stabParams.d_min;
    sysID(iOrder).ident(:, tooSmall) = [];
end

% Determine model orders range
n_min        = sysID(1).order;
n_step       = sysID(2).order - sysID(1).order;
n_max        = sysID(end).order;
modelOrders  = n_min : n_step : n_max;

% Prepare outputs
sysIDReduced = struct('order', {}, 'ident', {});
modeStab     = struct('fn', {}, 'dr', {}, 'mac', {}, 'globalMAC', {});

% Collect stable & unstable poles (for plotting)
stablePoles   = [];  % [order, freq, damp]
unstablePoles = [];  % [order, freq, damp]

%% ============= 3. Sort Each Order's Poles by Frequency =============
freqCell = cell(numel(sysID),1);
dampCell = cell(numel(sysID),1);
modeCell = cell(numel(sysID),1);

for iOrder = 1:numel(sysID)
    [sortedFreq, idx] = sort(sysID(iOrder).ident(1,:));
    sysID(iOrder).ident = sysID(iOrder).ident(:, idx);

    freqCell{iOrder} = sortedFreq;
    dampCell{iOrder} = sysID(iOrder).ident(2,:);
    modeCell{iOrder} = sysID(iOrder).ident(3:end,:);
end

%% ============= 4. Check Mode Stability Across Orders =============
for iOrder = 2:numel(sysID)
    [fStab, dStab, macStab] = compareModesBetweenOrders(...
        freqCell{iOrder},   freqCell{iOrder-1}, ...
        dampCell{iOrder},   dampCell{iOrder-1}, ...
        modeCell{iOrder},   modeCell{iOrder-1}, ...
        stabParams);

    % "Global" check: repeated detection across all orders
    globalMAC = false(1, size(sysID(iOrder).ident, 2));
    for iMode = 1:size(sysID(iOrder).ident, 2)
        thisFreq = sysID(iOrder).ident(1, iMode);
        countCloseFreq = 0;
        for jOrder = 1:numel(sysID)
            if jOrder == iOrder, continue; end
            freqDiff = abs(thisFreq - sysID(jOrder).ident(1,:));
            countCloseFreq = countCloseFreq + sum(freqDiff < stabParams.epsilon);
        end
        if countCloseFreq >= stabParams.NumMAC
            globalMAC(iMode) = true;
        end
    end

    finalCheck  = fStab & dStab & macStab & globalMAC;
    stableIdx   = finalCheck;
    unstableIdx = ~finalCheck;

    % Record stable poles
    stablePoles = [stablePoles; ...
        repmat(sysID(iOrder).order, sum(stableIdx),1), ...
        freqCell{iOrder}(stableIdx)', ...
        dampCell{iOrder}(stableIdx)'];

    % Record unstable poles
    unstablePoles = [unstablePoles; ...
        repmat(sysID(iOrder).order, sum(unstableIdx),1), ...
        freqCell{iOrder}(unstableIdx)', ...
        dampCell{iOrder}(unstableIdx)'];

    % Populate sysIDReduced
    sysIDReduced(iOrder).order = sysID(iOrder).order;
    sysIDReduced(iOrder).ident = sysID(iOrder).ident(:, stableIdx);

    % Store checks
    modeStab(iOrder).fn        = fStab;
    modeStab(iOrder).dr        = dStab;
    modeStab(iOrder).mac       = macStab;
    modeStab(iOrder).globalMAC = globalMAC;
end

%% ============= 5. Plot Depending on 'plotOption' + Save (if enabled) =============
% mainFig and mainAx are declared here so nested GUI callbacks can access
% them through the shared parent workspace (nested-function closure).
mainFig = [];   % Handle to the figure that hosts the interactive scatter
mainAx  = [];   % Handle to the left y-axis of that figure

switch lower(plotOption)
    case 'stab'
        % % ---------- FIGURE 1: Stabilisation Diagram (No Damping Color) ----------
        % fig1 = figure('Name','Stabilisation Diagram','NumberTitle','off',...
        %     'Color','white','Position',[300,300,1200,700]);
        %
        % % PSD on right y-axis
        % yyaxis right
        % hold on
        %
        % % ------- Flexible legend handling for PSD channels --------
        % if iscell(psdData.chan)
        %     legendStr = psdData.chan;
        %     nPSD = numel(psdData.chan);
        %     chVec = 1:nPSD;
        % else
        %     chVec = psdData.chan(:).';
        %     nPSD = numel(chVec);
        %     legendStr = arrayfun(@(ch) sprintf('PSD Ch. #%d', ch), chVec, 'UniformOutput', false);
        % end
        %
        % psdHandles = gobjects(nPSD,1);
        % for iCh = 1:nPSD
        %     psdHandles(iCh) = plot(psdData.fr_axis, psdData.data(:,chVec(iCh)), 'LineWidth',1);
        % end
        % set(gca,'YScale','log');
        % ylabel('Magnitude [dB]','FontSize',axisFontSize);
        % hold off
        %
        % % Poles on left y-axis
        % yyaxis left
        % hold on
        % plot(stablePoles(:,2), stablePoles(:,1), 'ro','MarkerSize',5);
        % plot(unstablePoles(:,2), unstablePoles(:,1), '.','Color',[0.2 0.2 0.2]);
        % hold off
        %
        % % Title, axes, legend
        % xlabel('Frequency [Hz]','FontSize',axisFontSize);
        % ylabel('Model Order','FontSize',axisFontSize);
        % title(plotTitle,'FontSize',mainFontSize,'FontWeight','bold');
        % grid on
        % set(gca,'FontSize',tickFontSize,'Box','off','YMinorTick','on','XMinorTick','on');
        % xlim([stabParams.f_min, stabParams.f_max]);
        % ylim([modelOrders(1)-5, modelOrders(end)+5]);
        %
        % % Create legend for PSD lines
        % lg = legend(psdHandles, legendStr, 'Location','north',"Orientation","Horizontal","FontSize",legendFontSize);

        % ---------- FIGURE 2: Stabilisation Diagram (Color-coded Damping) ----------
        fig2 = figure('Name','Stabilisation Diagram (Damping)','NumberTitle','off',...
            'Color','white','Position',[300,300,1200,700]);

        % PSD on right y-axis
        yyaxis right
        hold on
        % ------- Flexible legend handling for PSD channels --------
        if iscell(psdData.chan)
            legendStr2 = psdData.chan;
            nPSD2 = numel(psdData.chan);
            chVec2 = 1:nPSD2;
        else
            chVec2 = psdData.chan(:).';
            nPSD2 = numel(chVec2);
            legendStr2 = arrayfun(@(ch) sprintf('PSD Ch. #%d', ch), chVec2, 'UniformOutput', false);
        end

        psdHandles2 = gobjects(nPSD2,1);
        for iCh = 1:nPSD2
            psdHandles2(iCh) = plot(psdData.fr_axis, psdData.data(:,chVec2(iCh)), 'LineWidth',1);
        end
        set(gca,'YScale','log');
        ylabel('Magnitude [dB]','FontSize',axisFontSize);
        hold off

        % Poles on left y-axis (color by damping)
        yyaxis left
        hold on
        scatter(stablePoles(:,2), stablePoles(:,1), [], stablePoles(:,3), 'filled');
        caxis([stabParams.cbar_min, stabParams.cbar_max]);  % Damping colorbar range
        cb = colorbar;
        % Interpreter applied globally by applyTextInterpreter after plotting
        cb.Label.String          = 'Damping';
        cb.Label.FontSize        = axisFontSize;  % same as axis label
        cb.FontSize              = tickFontSize;  % colorbar numeric labels

        scatter(unstablePoles(:,2), unstablePoles(:,1), [], unstablePoles(:,3));
        hold off

        % Title, axes, legend
        xlabel('Frequency [Hz]','FontSize',axisFontSize);
        ylabel('Model Order','FontSize',axisFontSize);
        title(plotTitle,'FontSize',mainFontSize,'FontWeight','bold');
        grid on
        set(gca,'FontSize',tickFontSize,'Box','off','YMinorTick','on','XMinorTick','on');
        xlim([stabParams.f_min, stabParams.f_max]);
        ylim([modelOrders(1)-5, modelOrders(end)+5]);

        % Legend for PSD lines
        lg2 = legend(psdHandles2, legendStr2,'Location','north',...
            'Orientation','Horizontal','FontSize',legendFontSize); %#ok<NASGU>

        % Apply chosen text interpreter to every text element in fig2
        applyTextInterpreter(fig2, textInterpreter);

        % Store left-axis handle for the GUI callbacks (see Section 6 below)
        mainFig = fig2;
        yyaxis(fig2.CurrentAxes, 'left');
        mainAx  = fig2.CurrentAxes;

        % ----------------- Save if saveFigures = true -----------------
        if saveFigures
            saveFolder = 'ID_RESULTS';
            if ~exist(saveFolder,'dir'), mkdir(saveFolder); end

            if ~isempty(fileHandle)
                figNameBase = baseFileName + "_" + fileHandle;
            else
                figNameBase = baseFileName;
            end

            % Figure 1: No Damping Color  (fig1 is currently commented-out above)
            % exportgraphics(fig1, ...); etc.

            % Figure 2: Damping Color
            exportgraphics(fig2, saveFolder + "\" + figNameBase + '_stab_damping.pdf', ...
                'ContentType','vector','Resolution',600,'BackgroundColor','white');
            exportgraphics(fig2, saveFolder + "\" + figNameBase + '_stab_damping.png', ...
                'Resolution',300,'BackgroundColor','white');
            exportgraphics(fig2, saveFolder + "\" + figNameBase + '_stab_damping.tiff', ...
                'Resolution',600,'BackgroundColor','white');
            savefig(fig2, saveFolder + "\" + figNameBase + '_stab_damping.fig');
        end

    case 'clust'
        % ---------- FIGURE 3: Stabilisation Diagram (Cluster) ----------
        fig3 = figure('Name','Stabilisation Diagram (Cluster)','NumberTitle','off',...
            'Color','white','Position',[300,300,1000,600]);

        % PSD on right y-axis
        yyaxis right
        hold on
        if iscell(psdData.chan)
            legendStr3 = psdData.chan;
            nPSD3 = numel(psdData.chan);
            chVec3 = 1:nPSD3;
        else
            chVec3 = psdData.chan(:).';
            nPSD3 = numel(chVec3);
            legendStr3 = arrayfun(@(ch) sprintf('PSD Ch. #%d', ch), chVec3, 'UniformOutput', false);
        end
        psdHandles3 = gobjects(nPSD3,1);
        for iCh = 1:nPSD3
            psdHandles3(iCh) = plot(psdData.fr_axis, psdData.data(:,chVec3(iCh)), 'LineWidth',1);
        end
        set(gca,'YScale','log');
        ylabel('Magnitude [dB]','FontSize',axisFontSize);
        hold off

        % Poles on left y-axis
        yyaxis left
        hold on
        plot(stablePoles(:,2), stablePoles(:,1), 'ro','MarkerSize',5);
        plot(unstablePoles(:,2), unstablePoles(:,1), '.','Color',[0.2 0.2 0.2]);
        hold off

        xlabel('Frequency [Hz]','FontSize',axisFontSize);
        ylabel('Model Order','FontSize',axisFontSize);
        title(plotTitle,'FontSize',mainFontSize,'FontWeight','bold');
        grid on
        set(gca,'FontSize',tickFontSize,'Box','off','YMinorTick','on','XMinorTick','on');
        xlim([stabParams.f_min, stabParams.f_max]);
        ylim([modelOrders(1)-10, modelOrders(end)+0.5]);

        lg3 = legend(psdHandles3, legendStr3, 'Location','north',... %#ok<NASGU>
            'Orientation','Horizontal','FontSize',legendFontSize);
        applyTextInterpreter(fig3, textInterpreter);

        % ---------- FIGURE 4: Frequency-Damping Plot ----------
        fig4 = figure('Name','Frequency-Damping Plot','NumberTitle','off',...
            'Color','white','Position',[300,300,1000,600]); %#ok<NASGU>

        hold on
        plot(stablePoles(:,2), stablePoles(:,3), 'ro','MarkerSize',5);
        plot(unstablePoles(:,2), unstablePoles(:,3), '.','Color',[0.2 0.2 0.2]);
        hold off
        xlabel('Frequency [Hz]','FontSize',axisFontSize);
        ylabel('Damping','FontSize',axisFontSize);
        title([plotTitle ' (Freq-Damp)'],'FontSize',mainFontSize,'FontWeight','bold');
        grid on
        set(gca,'FontSize',tickFontSize,'Box','on','YMinorTick','on','XMinorTick','on');
        xlim([stabParams.f_min, stabParams.f_max]);
        ylim([stabParams.d_min, stabParams.d_max]);

        applyTextInterpreter(fig4, textInterpreter);

        % Store left-axis handle of fig3 (the order vs. freq plot) for GUI callbacks
        mainFig = fig3;
        yyaxis(fig3.CurrentAxes, 'left');
        mainAx  = fig3.CurrentAxes;

    case 'off'
        % No figures, no saving. selectedModes fields stay [].
end

%% ============= 6. Interactive Mode Selection (GUI) =============
% This section is skipped entirely when plotOption == 'off', or when
% interactive == false (tutorial repository: sysIDReduced(end) is taken
% as the highest-order selection instead -- see help text above).
% All GUI state is managed via nested functions below, which share the
% parent workspace variables through MATLAB's nested-function closure
% mechanism (no guidata / appdata required).

if interactive && ~strcmpi(plotOption, 'off') && ~isempty(mainFig) && ishandle(mainFig)

    % ---- Shared state for callbacks ----
    % selectedModes    : struct with .id [(2+nDOF)xN] and .order [1xN].
    % highlightHandles : 1xN array of graphics handles for the green markers.
    % Both are already declared in the outer workspace (selectedModes above,
    % highlightHandles below), so all nested functions share them directly.
    highlightHandles = gobjects(0, 1);  % Will grow as modes are selected

    % ---- Attach click handler to the main figure ----
    % Using WindowButtonDownFcn on the figure (rather than ButtonDownFcn on
    % a specific scatter object) provides a single, reliable entry point
    % regardless of which child graphics object intercepts the click.
    set(mainFig, 'WindowButtonDownFcn', @onFigureClicked);

    % ---- Build the control panel figure ----
    % Positioned to the left of the main figure to avoid overlap.
    ctrlFig = figure( ...
        'Name',        'Mode Selection Panel', ...
        'NumberTitle', 'off', ...
        'Color',       [0.94 0.94 0.94], ...
        'Position',    [50, 350, 360, 420], ...
        'MenuBar',     'none', ...
        'ToolBar',     'none', ...
        'Resize',      'off');

    % Instruction label
    uicontrol(ctrlFig, ...
        'Style',               'text', ...
        'String',              'Click stable poles on the main plot to select modes:', ...
        'Position',            [10 385 340 28], ...
        'HorizontalAlignment', 'left', ...
        'FontSize',            10, ...
        'BackgroundColor',     [0.94 0.94 0.94]);

    % Listbox: displays currently selected modes, one entry per row.
    % Format: "Freq: XX.XX Hz | Damp: X.XXX% | Ord: NNN"
    lbModes = uicontrol(ctrlFig, ...
        'Style',    'listbox', ...
        'Position', [10 110 340 268], ...
        'String',   {}, ...
        'FontName', 'Courier New', ...
        'FontSize', 10, ...
        'Value',    1);

    % "Delete Selected" button: removes the highlighted entry from the list
    % and the corresponding marker from the main plot.
    uicontrol(ctrlFig, ...
        'Style',    'pushbutton', ...
        'String',   'Delete Selected', ...
        'Position', [10 65 160 38], ...
        'FontSize', 10, ...
        'Callback', @onDeleteSelected);

    % "Done" button: resumes script execution and closes the control panel.
    uicontrol(ctrlFig, ...
        'Style',    'pushbutton', ...
        'String',   'Done', ...
        'Position', [185 65 165 38], ...
        'FontSize', 10, ...
        'FontWeight', 'bold', ...
        'Callback', @onDone);

    % Status label shown below the buttons
    uicontrol(ctrlFig, ...
        'Style',               'text', ...
        'String',              'Click on the stabilisation diagram, then press Done.', ...
        'Position',            [10 10 340 48], ...
        'HorizontalAlignment', 'center', ...
        'FontSize',            9, ...
        'ForegroundColor',     [0.3 0.3 0.3], ...
        'BackgroundColor',     [0.94 0.94 0.94]);

    % ---- Block execution until the user clicks Done (or closes ctrlFig) ----
    % uiwait suspends the calling function here.  It returns when uiresume
    % is called inside onDone, or when ctrlFig is closed via the window
    % manager (in which case selectedModes contains whatever was picked so far).
    uiwait(ctrlFig);

    % Remove the click handler from the main figure to restore normal behaviour
    if ishandle(mainFig)
        set(mainFig, 'WindowButtonDownFcn', '');
    end
end

%% ============= 7. Output: Cluster Data =============
clusterData = stablePoles;


% =========================================================================
% NESTED FUNCTIONS
% These are defined inside the main function body so that they share its
% workspace variables (stablePoles, selectedModes, highlightHandles,
% mainFig, mainAx, lbModes, ctrlFig) without any data-passing overhead.
% =========================================================================

    % -----------------------------------------------------------------
    % onFigureClicked  –  WindowButtonDownFcn for the main figure.
    % Fired on every mouse press over the figure.  Identifies the nearest
    % stable pole to the click position (in normalised axes coordinates)
    % and, if within the proximity threshold, adds it to the selection.
    % -----------------------------------------------------------------
    function onFigureClicked(src, ~)

        % Only respond to primary (left) button clicks
        if ~strcmp(get(src, 'SelectionType'), 'normal')
            return
        end

        % Switch the yyaxis context to 'left' before reading CurrentPoint
        % so that the returned coordinates are in [freq, order] space.
        yyaxis(mainAx, 'left');
        cp = get(mainAx, 'CurrentPoint');  % 2x3: top & bottom of view volume
        clickFreq  = cp(1, 1);             % x = frequency
        clickOrder = cp(1, 2);             % y = model order

        % Reject clicks outside the axes data limits
        xLim = get(mainAx, 'XLim');
        yLim = get(mainAx, 'YLim');
        if clickFreq  < xLim(1) || clickFreq  > xLim(2) || ...
           clickOrder < yLim(1) || clickOrder > yLim(2)
            return
        end

        % Guard: nothing to select if stablePoles is empty
        if isempty(stablePoles)
            return
        end

        % ---- Nearest-pole search in normalised axis space ----
        % Normalisation prevents frequency scale (e.g. 0–200 Hz) from
        % dominating over model-order scale (e.g. 2–60).
        xRange = xLim(2) - xLim(1);
        yRange = yLim(2) - yLim(1);

        dFreq  = (stablePoles(:,2) - clickFreq)  / xRange;
        dOrder = (stablePoles(:,1) - clickOrder) / yRange;
        dist   = sqrt(dFreq.^2 + dOrder.^2);

        [minDist, minIdx] = min(dist);

        % Proximity threshold: 3% of the normalised axes range.
        % Increase this value if selection feels unresponsive.
        if minDist > 0.03
            return
        end

        clickedPole = stablePoles(minIdx, :);  % [order, freq, damp]

        % ---- Duplicate check (freq + order uniqueness) ----
        if ~isempty(selMat)
            duplicate = any( ...
                abs(selMat(:,2) - clickedPole(2)) < 1e-9 & ...
                abs(selMat(:,1) - clickedPole(1)) < 1e-9 );
            if duplicate
                return
            end
        end

        % ---- Retrieve full ident column from sysIDReduced ----
        % sysIDReduced(k).ident contains all stable poles for order k.
        % Match by (order value, frequency) to extract the single column
        % corresponding to the clicked mode, including its mode shape.
        full_ident = [];  % fallback: empty if lookup fails
        for kk = 1:numel(sysIDReduced)
            if ~isempty(sysIDReduced(kk).order) && ...
               sysIDReduced(kk).order == clickedPole(1) && ...
               ~isempty(sysIDReduced(kk).ident)
                freq_row = sysIDReduced(kk).ident(1, :);   % Hz
                [~, col_match] = min(abs(freq_row - clickedPole(2)));
                full_ident = sysIDReduced(kk).ident(:, col_match);
                break
            end
        end
        if isempty(full_ident)
            % Fallback: assemble a minimal [freq; damp] column if the
            % lookup fails (should not happen in normal use).
            full_ident = clickedPole([2,3])';
        end

        % ---- Add mode to GUI state (selMat) and output matrix ----
        selMat(end+1, :) = clickedPole;          % [order, freq_Hz, damp]
        selectedModes.id    = [selectedModes.id,    full_ident];   % append column
        selectedModes.order = [selectedModes.order, clickedPole(1)]; % model order

        % Switch to the left y-axis of the main figure before plotting
        figure(mainFig);
        yyaxis(mainAx, 'left');
        hold(mainAx, 'on');
        hHL = plot(mainAx, clickedPole(2), clickedPole(1), 'o', ...
            'MarkerSize',      14, ...
            'MarkerEdgeColor', [0.0 0.75 0.0], ...  % bright green
            'MarkerFaceColor', 'none', ...
            'LineWidth',       2.0, ...
            'HandleVisibility','off');  % suppress legend entry
        hold(mainAx, 'off');

        highlightHandles(end+1) = hHL;

        % ---- Sort and refresh the control panel listbox ----
        sortAndUpdateDisplay();
    end % onFigureClicked


    % -----------------------------------------------------------------
    % onDeleteSelected  –  Callback for the "Delete Selected" button.
    % Removes the entry highlighted in the listbox, the corresponding
    % row from selectedModes, and the marker from the main plot.
    % -----------------------------------------------------------------
    function onDeleteSelected(~, ~)

        if isempty(selMat)
            return
        end

        selIdx = get(lbModes, 'Value');
        nModes = size(selMat, 1);

        if selIdx < 1 || selIdx > nModes
            return
        end

        % Remove the green marker from the main plot
        if ishandle(highlightHandles(selIdx))
            delete(highlightHandles(selIdx));
        end
        highlightHandles(selIdx) = [];

        % Remove from both state containers
        selMat(selIdx, :)        = [];   % remove row from GUI matrix
        selectedModes.id(:, selIdx) = [];   % remove column
        selectedModes.order(selIdx)  = [];   % remove order entry

        % Refresh the listbox; data are already sorted so no re-sort needed
        updateListbox();

        % Keep listbox selection within valid bounds after deletion
        newN = size(selMat, 1);
        if newN > 0
            set(lbModes, 'Value', min(selIdx, newN));
        else
            set(lbModes, 'Value', 1);
        end
    end % onDeleteSelected


    % -----------------------------------------------------------------
    % onDone  –  Callback for the "Done" button.
    % Resumes script execution (releases uiwait) and closes the control
    % panel.  selectedModes is already populated in the shared workspace.
    % -----------------------------------------------------------------
    function onDone(~, ~)
        % uiresume releases the uiwait call in the main function body.
        % The control panel figure is then closed; its WindowCloseRequestFcn
        % would also call uiresume, but calling it explicitly first is safer.
        if ishandle(ctrlFig)
            uiresume(ctrlFig);
            delete(ctrlFig);
        end
    end % onDone


    % -----------------------------------------------------------------
    % sortAndUpdateDisplay  –  Sorts selectedModes by frequency (ascending)
    % and synchronises highlightHandles to the same order, then refreshes
    % the listbox.  Called after every new selection.
    % -----------------------------------------------------------------
    function sortAndUpdateDisplay()

        if isempty(selMat)
            updateListbox();
            return
        end

        % Sort by frequency (column 2 of selMat) and apply the same
        % permutation to highlightHandles and selectedModes columns.
        [~, sortIdx]     = sort(selMat(:, 2));
        selMat           = selMat(sortIdx, :);
        highlightHandles = highlightHandles(sortIdx);
        selectedModes.id    = selectedModes.id(:, sortIdx);
        selectedModes.order = selectedModes.order(sortIdx);

        updateListbox();
    end % sortAndUpdateDisplay


    % -----------------------------------------------------------------
    % updateListbox  –  Rebuilds the listbox String cell array from the
    % current selectedModes matrix and pushes it to the uicontrol.
    % Kept separate from sortAndUpdateDisplay so that onDeleteSelected can
    % call it without triggering an unnecessary re-sort.
    % -----------------------------------------------------------------
    function updateListbox()
        % Display data comes from selMat (cheap matrix ops), not the struct.
        nModes = size(selMat, 1);

        if nModes == 0
            set(lbModes, 'String', {'(no modes selected)'}, 'Value', 1);
            return
        end

        entries = cell(nModes, 1);
        for k = 1:nModes
            % Format: "Freq: XXX.XX Hz | Damp: X.XXX% | Ord:  NN"
            entries{k} = sprintf('Freq: %7.3f Hz | Damp: %5.3f%% | Ord: %3d', ...
                selMat(k, 2), ...
                selMat(k, 3) * 100, ...
                selMat(k, 1));
        end

        % Guard Value against going out of range when items are added
        curVal = get(lbModes, 'Value');
        set(lbModes, 'String', entries, 'Value', min(curVal, nModes));
    end % updateListbox

end % === End of Main Function ===


%% ========================================================================
function [freqStable, dampStable, macStable] = ...
    compareModesBetweenOrders(freqHigh, freqLow, dampHigh, dampLow, ...
                              shapeHigh, shapeLow, stabParams)
% COMPAREMODESBETWEENORDERS
%   Checks frequency, damping, and MAC stability between consecutive model
%   orders. freqHigh/dampHigh correspond to the higher model order,
%   freqLow/dampLow to the lower model order.
%
%   Returns logical arrays (freqStable, dampStable, macStable) for each mode
%   in freqHigh/dampHigh.

freqStable = false(size(freqHigh));
dampStable = false(size(dampHigh));
macStable  = false(size(freqHigh));

% If the lower model order is empty, do nothing
if isempty(freqLow)
    return
end

for iMode = 1:length(freqHigh)
    % Frequency stability (relative difference)
    relFreqDiff  = abs(freqLow - freqHigh(iMode)) / freqHigh(iMode);
    freqStable(iMode) = min(relFreqDiff) < stabParams.dfr;

    % Damping stability (relative difference)
    relDampDiff  = abs(dampLow - dampHigh(iMode)) / dampHigh(iMode);
    dampStable(iMode) = min(relDampDiff) < stabParams.dz;

    % MAC check
    macVals = compute_mac(shapeHigh(:,iMode), shapeLow);
    if max(macVals) > stabParams.dMAC
        macStable(iMode) = true;
    end
end

end


%% ========================================================================
function macVal = compute_mac(vec, mat)
% COMPUTE_MAC
%   Computes the Modal Assurance Criterion (MAC) between a single vector 'vec'
%   and each column of 'mat'. MAC is defined by:
%
%       MAC = (|v' * u|^2) / ((v' * v) * (u' * u))
%
%   Returns a row of MAC values for each column in 'mat'.

if size(vec,2) > 1
    vec = vec(:);
end

macVal = zeros(1, size(mat,2));
for j = 1:size(mat,2)
    u = mat(:,j);
    numerator   = abs(vec' * u)^2;
    denominator = (vec' * vec) * (u' * u);
    if denominator > 0
        macVal(j) = numerator / denominator;
    else
        macVal(j) = 0;
    end
end

end


%% ========================================================================
function applyTextInterpreter(fig, interp)
% APPLYTEXTINTERPRETER  Apply a text interpreter to every text element of a figure.
%
%   Covers: titles, axis labels, legend entries, annotation text objects
%   (via the 'Interpreter' property), plus axes tick labels and colorbar
%   tick labels (via 'TickLabelInterpreter').
%
%   Called once per figure after all plotting is complete so that elements
%   added later (e.g., colorbars, legends) are also captured.

% All graphics objects that carry an 'Interpreter' property
set(findall(fig, '-property', 'Interpreter'), 'Interpreter', interp);

% Axes tick labels (separate property from text Interpreter)
set(findall(fig, 'Type', 'axes'), 'TickLabelInterpreter', interp);

% Colorbar tick labels
cbObjs = findall(fig, 'Type', 'colorbar');
for iCB = 1:numel(cbObjs)
    cbObjs(iCB).TickLabelInterpreter = interp;
end

end