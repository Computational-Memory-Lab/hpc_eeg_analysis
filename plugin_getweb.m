function [pluginTmp, eeglabVersionStatus] = plugin_getweb(varargin) %#ok<INUSD>
% PLUGIN_GETWEB - Offline-safe fallback for EEGLAB startup checks.
%
% In cluster batch jobs, intermittent network/filesystem issues can make
% EEGLAB's remote plugin metadata fetch fail and crash startup. Returning a
% minimal placeholder struct keeps startup robust in non-interactive runs.

pluginTmp = struct('name', 'offline_placeholder', 'version', '0');
eeglabVersionStatus = struct();
end
