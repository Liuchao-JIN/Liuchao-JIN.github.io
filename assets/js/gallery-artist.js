(function () {
  "use strict";

  var root = document.querySelector("[data-artist-gallery]");
  if (!root) return;

  var cards = Array.prototype.slice.call(root.querySelectorAll("[data-gallery-theme]"));
  var filters = Array.prototype.slice.call(root.querySelectorAll("[data-gallery-filter]"));
  var lightbox = document.querySelector("[data-gallery-lightbox]");
  var lightboxImage = document.querySelector("[data-gallery-lightbox-image]");
  var lightboxCaption = document.querySelector("[data-gallery-lightbox-caption]");
  var lightboxMeta = document.querySelector("[data-gallery-lightbox-meta]");
  var closeButton = document.querySelector("[data-gallery-close]");
  var previousButton = document.querySelector("[data-gallery-prev]");
  var nextButton = document.querySelector("[data-gallery-next]");
  var activeItems = [];
  var activeIndex = -1;
  var lastTrigger = null;
  var touchStartX = null;

  function refreshActiveItems() {
    activeItems = cards
      .filter(function (card) { return !card.hidden; })
      .map(function (card) { return card.querySelector("[data-gallery-item]"); });
  }

  function renderItem(index) {
    if (!activeItems.length) return;
    activeIndex = (index + activeItems.length) % activeItems.length;
    var item = activeItems[activeIndex];
    lightboxImage.src = item.getAttribute("data-src");
    lightboxImage.alt = item.getAttribute("data-caption") || "Gallery photograph";
    lightboxCaption.textContent = item.getAttribute("data-caption") || "";
    lightboxMeta.textContent = item.getAttribute("data-meta") || "";
  }

  function openLightbox(item) {
    refreshActiveItems();
    activeIndex = activeItems.indexOf(item);
    if (activeIndex < 0) return;
    lastTrigger = item;
    renderItem(activeIndex);
    lightbox.classList.add("is-open");
    lightbox.setAttribute("aria-hidden", "false");
    document.body.classList.add("gallery-lightbox-open");
    closeButton.focus();
  }

  function closeLightbox() {
    lightbox.classList.remove("is-open");
    lightbox.setAttribute("aria-hidden", "true");
    document.body.classList.remove("gallery-lightbox-open");
    lightboxImage.removeAttribute("src");
    if (lastTrigger) lastTrigger.focus();
  }

  filters.forEach(function (button) {
    button.addEventListener("click", function () {
      var filter = button.getAttribute("data-gallery-filter");
      filters.forEach(function (candidate) {
        var selected = candidate === button;
        candidate.classList.toggle("is-active", selected);
        candidate.setAttribute("aria-pressed", selected ? "true" : "false");
      });
      cards.forEach(function (card) {
        card.hidden = filter !== "all" && card.getAttribute("data-gallery-theme") !== filter;
      });
      refreshActiveItems();
    });
  });

  cards.forEach(function (card) {
    var item = card.querySelector("[data-gallery-item]");
    item.addEventListener("click", function () { openLightbox(item); });
  });

  closeButton.addEventListener("click", closeLightbox);
  previousButton.addEventListener("click", function () { renderItem(activeIndex - 1); });
  nextButton.addEventListener("click", function () { renderItem(activeIndex + 1); });

  lightbox.addEventListener("click", function (event) {
    if (event.target === lightbox) closeLightbox();
  });

  lightbox.addEventListener("touchstart", function (event) {
    touchStartX = event.changedTouches[0].clientX;
  }, { passive: true });

  lightbox.addEventListener("touchend", function (event) {
    if (touchStartX === null) return;
    var distance = event.changedTouches[0].clientX - touchStartX;
    if (Math.abs(distance) > 55) renderItem(activeIndex + (distance < 0 ? 1 : -1));
    touchStartX = null;
  }, { passive: true });

  document.addEventListener("keydown", function (event) {
    if (!lightbox.classList.contains("is-open")) return;
    if (event.key === "Escape") closeLightbox();
    if (event.key === "ArrowLeft") renderItem(activeIndex - 1);
    if (event.key === "ArrowRight") renderItem(activeIndex + 1);
  });

  refreshActiveItems();
})();
