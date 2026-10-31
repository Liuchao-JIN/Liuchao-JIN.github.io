(function () {
  "use strict";

  function onReady(callback) {
    if (document.readyState === "loading") {
      document.addEventListener("DOMContentLoaded", callback);
    } else {
      callback();
    }
  }

  function normalizeText(value) {
    return (value || "").replace(/\s+/g, " ").trim();
  }

  function normalizeLocalUrl(value) {
    return (value || "").replace(/^http:\/\/liuchao-jin\.github\.io/i, "");
  }

  function cloneNodesBetween(startNode, endNode) {
    var nodes = [];
    var current = startNode ? startNode.nextSibling : null;

    while (current && current !== endNode) {
      var isEmptyParagraph = current.nodeType === 1 &&
        current.tagName === "P" &&
        !normalizeText(current.textContent);

      if (!isEmptyParagraph) {
        nodes.push(current.cloneNode(true));
      }
      current = current.nextSibling;
    }
    return nodes;
  }

  function titleBefore(firstBreak, container) {
    var parts = [];
    var current = container.firstChild;

    while (current && current !== firstBreak) {
      parts.push(current.textContent || "");
      current = current.nextSibling;
    }
    return normalizeText(parts.join(" "));
  }

  function appendNodes(container, nodes) {
    nodes.forEach(function (node) {
      container.appendChild(node);
    });
  }

  function createAuthors(authorNodes) {
    var paragraph = document.createElement("p");
    paragraph.className = "publication-authors";
    appendNodes(paragraph, authorNodes);

    var authorText = normalizeText(paragraph.textContent);
    var authorCount = authorText ? authorText.split(",").length : 0;

    if (authorText.length <= 285 && authorCount <= 11) {
      return paragraph;
    }

    var firstAuthor = normalizeText(authorText.split(",")[0]);
    var details = document.createElement("details");
    var summary = document.createElement("summary");
    var fullList = document.createElement("p");

    details.className = "publication-authors publication-authors--long";
    summary.textContent = firstAuthor + ", …, Liuchao Jin, … · Show all " + authorCount + " authors";
    appendNodes(fullList, authorNodes.map(function (node) {
      return node.cloneNode(true);
    }));
    details.appendChild(summary);
    details.appendChild(fullList);
    return details;
  }

  function createAction(link, title, toast) {
    var originalLabel = normalizeText(link.textContent).toLowerCase();
    var href = normalizeLocalUrl(link.getAttribute("href"));

    if (originalLabel.indexOf("bibtex") !== -1) {
      var button = document.createElement("button");
      button.type = "button";
      button.textContent = "Copy BibTeX";
      button.setAttribute("aria-label", "Copy BibTeX for " + title);
      button.addEventListener("click", function () {
        fetch(href, { credentials: "same-origin" })
          .then(function (response) {
            if (!response.ok) {
              throw new Error("BibTeX request failed");
            }
            return response.text();
          })
          .then(copyText)
          .then(function () {
            showToast(toast, "BibTeX copied");
          })
          .catch(function () {
            window.open(href, "_blank", "noopener");
            showToast(toast, "Opened BibTeX file");
          });
      });
      return button;
    }

    var action = document.createElement("a");
    var isWeb = originalLabel === "web" || originalLabel === "doi";
    action.href = href;
    action.textContent = isWeb ? "DOI" : normalizeText(link.textContent);
    action.target = "_blank";
    action.rel = "noopener";
    action.setAttribute("aria-label", action.textContent + " for " + title);

    if (isWeb) {
      action.className = "publication-action--primary";
    }
    return action;
  }

  function copyText(text) {
    if (navigator.clipboard && window.isSecureContext) {
      return navigator.clipboard.writeText(text);
    }

    return new Promise(function (resolve, reject) {
      var input = document.createElement("textarea");
      input.value = text;
      input.setAttribute("readonly", "");
      input.style.position = "fixed";
      input.style.opacity = "0";
      document.body.appendChild(input);
      input.select();

      try {
        document.execCommand("copy");
        resolve();
      } catch (error) {
        reject(error);
      } finally {
        document.body.removeChild(input);
      }
    });
  }

  var toastTimer;
  function showToast(toast, message) {
    if (!toast) {
      return;
    }
    window.clearTimeout(toastTimer);
    toast.textContent = message;
    toast.classList.add("is-visible");
    toastTimer = window.setTimeout(function () {
      toast.classList.remove("is-visible");
    }, 1600);
  }

  function slugify(value) {
    return normalizeText(value)
      .toLowerCase()
      .replace(/[^a-z0-9]+/g, "-")
      .replace(/^-|-$/g, "")
      .slice(0, 72);
  }

  function enhanceTable(table, year, toast) {
    var list = document.createElement("div");
    var articles = [];
    var usedIds = {};

    list.className = "publications-list";
    list.dataset.publicationYear = year;

    Array.prototype.forEach.call(table.querySelectorAll("tbody > tr"), function (row) {
      var image = row.querySelector("td:first-child img");
      var information = row.querySelector("td:nth-child(2)");

      if (!image || !information) {
        return;
      }

      var breaks = information.querySelectorAll("br");
      var journalElement = information.querySelector("p em");
      var links = Array.prototype.slice.call(information.querySelectorAll("a"));

      if (!breaks.length || !journalElement) {
        return;
      }

      var title = titleBefore(breaks[0], information);
      var authorNodes = cloneNodesBetween(breaks[0], breaks[1]);
      var authors = createAuthors(authorNodes);
      var journal = normalizeText(journalElement.textContent);
      var article = document.createElement("article");
      var figure = document.createElement("div");
      var publicationImage = image.cloneNode(true);
      var meta = document.createElement("div");
      var yearBadge = document.createElement("span");
      var journalName = document.createElement("em");
      var heading = document.createElement("h3");
      var actions = document.createElement("div");
      var baseId = "publication-" + slugify(title);
      var articleId = baseId;
      var duplicateIndex = 2;

      while (usedIds[articleId]) {
        articleId = baseId + "-" + duplicateIndex;
        duplicateIndex += 1;
      }
      usedIds[articleId] = true;

      article.id = articleId;
      article.className = "publication-card";
      article.dataset.publicationYear = year;
      article.dataset.publicationSearch = normalizeText(
        title + " " + authors.textContent + " " + journal + " " + year
      ).toLowerCase();

      figure.className = "publication-figure";
      publicationImage.removeAttribute("width");
      publicationImage.removeAttribute("style");
      publicationImage.loading = "lazy";
      publicationImage.decoding = "async";
      if (!publicationImage.getAttribute("alt")) {
        publicationImage.alt = "Representative figure for " + title;
      }
      figure.appendChild(publicationImage);

      meta.className = "publication-meta";
      yearBadge.textContent = year;
      journalName.textContent = journal;
      meta.appendChild(yearBadge);
      meta.appendChild(journalName);

      heading.className = "publication-title";
      heading.textContent = title;

      actions.className = "publication-actions";
      actions.setAttribute("aria-label", "Links for " + title);
      links.forEach(function (link) {
        actions.appendChild(createAction(link, title, toast));
      });

      article.appendChild(figure);
      article.appendChild(meta);
      article.appendChild(heading);
      article.appendChild(authors);
      article.appendChild(actions);
      list.appendChild(article);
      articles.push(article);
    });

    table.replaceWith(list);
    return { list: list, articles: articles };
  }

  function updateUrl(year, query) {
    try {
      var url = new URL(window.location.href);
      if (year && year !== "all") {
        url.searchParams.set("year", year);
      } else {
        url.searchParams.delete("year");
      }
      if (query) {
        url.searchParams.set("q", query);
      } else {
        url.searchParams.delete("q");
      }
      window.history.replaceState({}, "", url.toString());
    } catch (error) {
      // Filtering still works when URL state cannot be updated.
    }
  }

  onReady(function () {
    var page = document.querySelector("[data-publications-page]");
    if (!page) {
      return;
    }

    var controls = page.querySelector("[data-publications-controls]");
    var search = page.querySelector("[data-publications-search]");
    var yearSelect = page.querySelector("[data-publications-year-select]");
    var yearButtons = Array.prototype.slice.call(
      page.querySelectorAll("[data-publications-year]")
    );
    var resultCount = page.querySelector("[data-publications-result-count]");
    var summary = page.querySelector("[data-publications-summary]");
    var toast = page.querySelector("[data-publications-toast]");
    var groups = [];
    var allArticles = [];

    Array.prototype.forEach.call(page.querySelectorAll("table.publications-list"), function (table) {
      var year = table.dataset.publicationYear;
      var heading = page.querySelector("#publications-" + year);
      var enhanced = enhanceTable(table, year, toast);
      var count = document.createElement("span");

      count.className = "publication-year-count";
      count.textContent = enhanced.articles.length + " publications";
      heading.appendChild(count);

      groups.push({
        year: year,
        heading: heading,
        list: enhanced.list,
        articles: enhanced.articles,
        count: count
      });
      allArticles = allArticles.concat(enhanced.articles);
    });

    var empty = document.createElement("p");
    empty.className = "publications-empty";
    empty.textContent = "No publications match the current filters.";
    empty.hidden = true;
    page.insertBefore(empty, toast);

    var params = new URLSearchParams(window.location.search);
    var availableYears = groups.map(function (group) { return group.year; });
    var activeYear = params.get("year") || "all";

    if (activeYear !== "all" && availableYears.indexOf(activeYear) === -1) {
      activeYear = "all";
    }
    search.value = params.get("q") || "";

    function renderFilters(writeUrl) {
      var query = normalizeText(search.value).toLowerCase();
      var visibleTotal = 0;

      groups.forEach(function (group) {
        var visibleInYear = 0;

        group.articles.forEach(function (article) {
          var matchesYear = activeYear === "all" || article.dataset.publicationYear === activeYear;
          var matchesSearch = !query || article.dataset.publicationSearch.indexOf(query) !== -1;
          var visible = matchesYear && matchesSearch;
          article.hidden = !visible;
          if (visible) {
            visibleInYear += 1;
            visibleTotal += 1;
          }
        });

        group.heading.hidden = visibleInYear === 0;
        group.list.hidden = visibleInYear === 0;
        group.count.textContent = visibleInYear === group.articles.length
          ? group.articles.length + " publications"
          : visibleInYear + " of " + group.articles.length + " publications";
      });

      yearButtons.forEach(function (button) {
        var selected = button.dataset.publicationsYear === activeYear;
        button.setAttribute("aria-pressed", selected ? "true" : "false");
      });
      yearSelect.value = activeYear;
      resultCount.textContent = visibleTotal + (visibleTotal === 1 ? " publication" : " publications");
      empty.hidden = visibleTotal !== 0;

      if (writeUrl) {
        updateUrl(activeYear, query);
      }
    }

    yearButtons.forEach(function (button) {
      button.addEventListener("click", function () {
        activeYear = button.dataset.publicationsYear;
        renderFilters(true);
      });
    });

    yearSelect.addEventListener("change", function () {
      activeYear = yearSelect.value;
      renderFilters(true);
    });

    search.addEventListener("input", function () {
      renderFilters(true);
    });

    summary.textContent = allArticles.length + " publications · " +
      availableYears[availableYears.length - 1] + "–" + availableYears[0];
    controls.hidden = false;
    page.classList.add("is-enhanced");
    renderFilters(false);
  });
})();
