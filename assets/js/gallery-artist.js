(function () {
  "use strict";

  var root = document.querySelector("[data-artist-gallery]");
  var lightbox = document.querySelector("[data-gallery-lightbox]");
  if (!root || !lightbox) return;

  var items = Array.prototype.slice.call(root.querySelectorAll("[data-gallery-item]"));
  if (!items.length) return;
  var lightboxImage = document.querySelector("[data-gallery-lightbox-image]");
  var lightboxCaption = document.querySelector("[data-gallery-lightbox-caption]");
  var lightboxMeta = document.querySelector("[data-gallery-lightbox-meta]");
  var closeButton = document.querySelector("[data-gallery-close]");
  var previousButton = document.querySelector("[data-gallery-prev]");
  var nextButton = document.querySelector("[data-gallery-next]");
  var activeIndex = -1;
  var lastTrigger = null;
  var touchStartX = null;

  function renderItem(index) {
    activeIndex = (index + items.length) % items.length;
    var item = items[activeIndex];
    lightboxImage.src = item.getAttribute("data-src");
    lightboxImage.alt = item.querySelector("img").alt;
    lightboxCaption.textContent = item.getAttribute("data-caption") || "";
    lightboxMeta.textContent = item.getAttribute("data-meta") || "";
  }

  function openLightbox(item) {
    lastTrigger = item;
    renderItem(items.indexOf(item));
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
    touchStartX = null;
    if (lastTrigger) lastTrigger.focus();
  }

  items.forEach(function (item) {
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
    if (event.key === "ArrowLeft" || event.key === "ArrowRight") {
      event.preventDefault();
      renderItem(activeIndex + (event.key === "ArrowLeft" ? -1 : 1));
    }
    if (event.key === "Tab" && event.shiftKey && document.activeElement === closeButton) {
      event.preventDefault();
      nextButton.focus();
    } else if (event.key === "Tab" && !event.shiftKey && document.activeElement === nextButton) {
      event.preventDefault();
      closeButton.focus();
    }
  });
})();
