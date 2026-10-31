(function () {
  "use strict";

  function onReady(callback) {
    if (document.readyState === "loading") {
      document.addEventListener("DOMContentLoaded", callback);
    } else {
      callback();
    }
  }

  onReady(function () {
    var page = document.querySelector("[data-cv-page]");
    var printButton = document.querySelector("[data-cv-print]");
    var originalDetailStates = [];

    if (!page || !printButton) {
      return;
    }

    printButton.addEventListener("click", function () {
      window.print();
    });

    window.addEventListener("beforeprint", function () {
      originalDetailStates = Array.prototype.map.call(page.querySelectorAll("details"), function (details) {
        var wasOpen = details.open;
        details.open = true;
        return { details: details, wasOpen: wasOpen };
      });
    });

    window.addEventListener("afterprint", function () {
      originalDetailStates.forEach(function (item) {
        item.details.open = item.wasOpen;
      });
      originalDetailStates = [];
    });
  });
})();
