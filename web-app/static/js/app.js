(function ($) {
    const state = {
        config: null,
        lastResults: null,
    };

    function setStatus(selector, message, kind) {
        const element = $(selector);
        element.text(message || "");
        element.removeClass("status-error status-success");
        if (kind === "error") element.addClass("status-error");
        if (kind === "success") element.addClass("status-success");
    }

    function renderLoadedDatabase() {
        const panel = $("#db-info-panel");
        panel.empty();

        if (!state.config || !state.config.db_file) {
            panel.append('<div class="small text-body-secondary">Start the server with --db-file /path/to/database.h5</div>');
            return;
        }

        const rows = [];
        if (state.config.fp_type) rows.push(["Fingerprint type", state.config.fp_type]);
        if (state.config.fp_bits_size != null) rows.push(["FP size (bits)", state.config.fp_bits_size]);
        if (state.config.num_mols != null) rows.push(["Total compounds", Number(state.config.num_mols).toLocaleString()]);
        if (state.config.cluster_threshold != null) rows.push(["Cluster threshold", state.config.cluster_threshold]);
        rows.push(["Search mode", state.config.search_mode || "memory"]);
        rows.push(["Results mode", state.config.results_mode || "ids_only"]);

        if (state.config.db_info && Object.keys(state.config.db_info).length) {
            Object.entries(state.config.db_info).forEach(([k, v]) => {
                rows.push([k, typeof v === "object" ? JSON.stringify(v) : String(v)]);
            });
        }

        const table = $('<table class="info-table"></table>');
        rows.forEach(([label, value]) => {
            table.append(`<tr><th>${label}</th><td>${value}</td></tr>`);
        });
        panel.append(table);
    }

    function csvEscape(value) {
        if (value == null) return "";
        const str = String(value);
        if (str.includes(",") || str.includes('"') || str.includes("\n")) {
            return '"' + str.replace(/"/g, '""') + '"';
        }
        return str;
    }

    function downloadCSV() {
        if (!state.lastResults) return;
        const hasSmiles = state.lastResults.some(r => r.hits && r.hits.some(h => h.smiles));
        const headers = ["Query", "Query SMILES", "Compound ID", "Score"];
        if (hasSmiles) headers.push("SMILES");
        const lines = [headers.join(",")];
        state.lastResults.forEach((result, qi) => {
            if (!result.hits) return;
            result.hits.forEach(hit => {
                const row = [
                    qi + 1,
                    csvEscape(result.query_smiles),
                    csvEscape(hit.compound_id),
                    hit.score.toFixed(4),
                ];
                if (hasSmiles) row.push(csvEscape(hit.smiles || ""));
                lines.push(row.join(","));
            });
        });
        const blob = new Blob([lines.join("\n")], { type: "text/csv;charset=utf-8;" });
        const url = URL.createObjectURL(blob);
        const a = document.createElement("a");
        a.href = url;
        a.download = "uffpsim_results.csv";
        a.click();
        URL.revokeObjectURL(url);
    }

    function renderResults(payload) {
        const resultsContainer = $("#results-container");
        resultsContainer.empty();
        $("#download-csv").hide();
        state.lastResults = null;

        const results = payload.results || [];
        const totalHits = results.reduce((sum, r) => sum + (r.hits ? r.hits.length : 0), 0);

        if (!results.length) {
            resultsContainer.append('<div class="empty-state">No results yet. Run a search to populate this area.</div>');
            $("#results-summary").text("");
            return;
        }

        $("#results-summary").text(`${results.length} query result(s) — ${totalHits} total hit(s)`);
        state.lastResults = results;

        const hasSmiles = results.some(r => r.hits && r.hits.some(h => h.smiles));
        const hasImage  = results.some(r => r.hits && r.hits.some(h => h.image));

        results.forEach((result, queryIndex) => {
            const queryCard = $('<article class="query-card"></article>');

            const queryImageHtml = result.query_image
                ? `<img class="molecule-image" src="${result.query_image}" alt="Query ${queryIndex + 1}" style="max-width:140px">`
                : "";

            queryCard.append(`
                <div class="query-card-header">
                    <div class="d-flex gap-3 align-items-center flex-wrap">
                        ${queryImageHtml}
                        <div class="flex-grow-1">
                            <h3 class="h5 mb-1">Query ${queryIndex + 1}</h3>
                            <div class="smiles-block">${result.query_smiles}</div>
                        </div>
                        <span class="badge text-bg-primary">${result.hits ? result.hits.length : 0} hit(s)</span>
                    </div>
                </div>
            `);

            let bodyHtml = "";
            if (result.error) {
                bodyHtml = `<div class="alert alert-warning mb-0">${result.error}</div>`;
            } else if (!result.hits || !result.hits.length) {
                bodyHtml = '<div class="empty-state">No hits matched the requested threshold.</div>';
            } else {
                let headerCols = "<th>#</th><th>Compound ID</th><th>Score</th>";
                if (hasSmiles) headerCols += "<th>SMILES</th>";
                if (hasImage)  headerCols += "<th>Structure</th>";

                const rowsHtml = result.hits.map((hit, i) => {
                    let cols = `<td>${i + 1}</td><td><code>${hit.compound_id}</code></td><td>${hit.score.toFixed(4)}</td>`;
                    if (hasSmiles) cols += `<td><span class="smiles-block">${hit.smiles || "—"}</span></td>`;
                    if (hasImage)  cols += hit.image
                        ? `<td><img class="hit-image" src="${hit.image}" alt="${hit.compound_id}"></td>`
                        : `<td>—</td>`;
                    return `<tr>${cols}</tr>`;
                }).join("");

                bodyHtml = `<div class="table-responsive">
                    <table class="results-table">
                        <thead><tr>${headerCols}</tr></thead>
                        <tbody>${rowsHtml}</tbody>
                    </table>
                </div>`;
            }

            queryCard.append(`<div class="query-card-body">${bodyHtml}</div>`);
            resultsContainer.append(queryCard);
        });

        if (totalHits > 0) $("#download-csv").show();
    }

    function loadConfig() {
        return $.getJSON("/api/config").done((payload) => {
            state.config = payload;
            renderLoadedDatabase();
        });
    }

    function submitSearch(event) {
        event.preventDefault();
        setStatus("#search-status", "Searching…", null);
        $("#results-summary").text("");

        if (!state.config || !state.config.db_file) {
            setStatus("#search-status", "Server was not started with a database file.", "error");
            return;
        }

        const payload = {
            threshold: Number($("#search-threshold").val()),
            limit_by: Number($("#search-limit").val()),
            smiles_text: $("#smiles-input").val(),
        };

        $.ajax({
            url: "/api/search",
            method: "POST",
            contentType: "application/json",
            data: JSON.stringify(payload),
        }).done((response) => {
            renderResults(response);
            setStatus("#search-status", `Search complete for ${response.results.length} query entries.`, "success");
        }).fail((xhr) => {
            const message = xhr.responseJSON && xhr.responseJSON.error ? xhr.responseJSON.error : "Search failed.";
            setStatus("#search-status", message, "error");
            renderResults({ results: [] });
        });
    }

    $(function () {
        loadConfig().then(() => {
            renderResults({ results: [] });
        });

        $("#search-form").on("submit", submitSearch);
        $("#download-csv").on("click", downloadCSV);
    });
})(jQuery);
