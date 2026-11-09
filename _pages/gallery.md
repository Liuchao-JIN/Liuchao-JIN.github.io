---
layout: gallery
title: "Gallery"
permalink: /gallery/
excerpt: "A visual field notebook of travel, academic journeys, architecture, landscape, and quiet observation."
author_profile: false
hide_title: true
---

<div class="artist-gallery" data-artist-gallery>
  <section class="artist-gallery__hero" aria-labelledby="gallery-title">
    <div class="artist-gallery__hero-copy">
      <p class="artist-gallery__eyebrow">A VISUAL FIELD NOTEBOOK · 2017—2026</p>
      <h1 id="gallery-title">Light, distance,<br>and the spaces between.</h1>
      <p class="artist-gallery__lede">A personal archive gathered between research journeys: mountains emerging from cloud, cities after dark, rooms made quiet by art, and brief encounters with unfamiliar places.</p>
      <a class="artist-gallery__scroll" href="#gallery-archive">Enter the archive <span aria-hidden="true">↓</span></a>
    </div>

    <div class="artist-gallery__poster" aria-label="Selected photographs">
      <figure class="artist-gallery__poster-main">
        <img src="{{ '/images/gallery/2025-06-07_Nepal_Himalayan_Valley_0804.jpg' | relative_url }}" alt="Himalayan ridges and valleys in Nepal" width="1800" height="1350" fetchpriority="high">
        <figcaption>Nepal · 2025</figcaption>
      </figure>
      <figure class="artist-gallery__poster-top">
        <img src="{{ '/images/gallery/2024-07-14_London_UK_Contemporary_Gallery_1513.jpg' | relative_url }}" alt="A quiet contemporary gallery in London" width="1800" height="1350" fetchpriority="high">
        <figcaption>London · 2024</figcaption>
      </figure>
      <figure class="artist-gallery__poster-bottom">
        <img src="{{ '/images/gallery/2025-12-22_Nevada_USA_Red_Rock_Portrait_1353.jpg' | relative_url }}" alt="Portrait among the red rocks of Nevada" width="1350" height="1800" fetchpriority="high">
        <figcaption>Nevada · 2025</figcaption>
      </figure>
      <span class="artist-gallery__poster-mark" aria-hidden="true">FIELD<br>NOTES</span>
    </div>
  </section>

  <section class="artist-gallery__statement" aria-label="Artist statement">
    <span class="artist-gallery__statement-number">01</span>
    <p>These are not records of arrival. They are studies of atmosphere—how a place holds light, how scale changes attention, and how travel turns ordinary moments into composition.</p>
  </section>

  <section class="artist-gallery__archive" id="gallery-archive" aria-labelledby="archive-title">
    <header class="artist-gallery__archive-head">
      <div>
        <p class="artist-gallery__eyebrow">THE ARCHIVE</p>
        <h2 id="archive-title">Collected observations</h2>
      </div>
      <p class="artist-gallery__archive-note">Select a photograph to enter the frame.</p>
    </header>

    <div class="artist-gallery__filters" role="group" aria-label="Filter photographs by theme">
      <button class="is-active" type="button" data-gallery-filter="all" aria-pressed="true">All <span>{{ site.data.gallery | size }}</span></button>
      <button type="button" data-gallery-filter="landscape" aria-pressed="false">Landscape</button>
      <button type="button" data-gallery-filter="city" aria-pressed="false">City</button>
      <button type="button" data-gallery-filter="culture" aria-pressed="false">Culture</button>
      <button type="button" data-gallery-filter="people" aria-pressed="false">People</button>
      <button type="button" data-gallery-filter="journey" aria-pressed="false">Journey</button>
    </div>

    <div class="artist-gallery__grid" data-gallery-grid aria-live="polite">
      {% for photo in site.data.gallery %}
      <figure class="artist-gallery__card artist-gallery__card--{{ photo.orientation }}" data-gallery-theme="{{ photo.theme }}">
        <button type="button"
                class="artist-gallery__open"
                data-gallery-item
                data-src="{{ photo.src | relative_url }}"
                data-caption="{{ photo.caption | escape }}"
                data-meta="{{ photo.place | escape }} · {{ photo.year }}"
                aria-label="Open photograph: {{ photo.caption | escape }}">
          <span class="artist-gallery__image-wrap">
            <img src="{{ photo.src | relative_url }}"
                 alt="{{ photo.alt | escape }}"
                 loading="lazy"
                 decoding="async">
            <span class="artist-gallery__image-index" aria-hidden="true">{{ photo.id }}</span>
          </span>
          <span class="artist-gallery__caption">
            <span>{{ photo.place }}</span>
            <span>{{ photo.year }}</span>
          </span>
        </button>
      </figure>
      {% endfor %}
    </div>
  </section>

  <section class="artist-gallery__closing" aria-label="Closing note">
    <p class="artist-gallery__eyebrow">A CONTINUING PRACTICE</p>
    <p class="artist-gallery__closing-title">The archive remains open.</p>
    <p>New places will enter slowly—one frame, one encounter, one remembered light at a time.</p>
    <a href="{{ '/amateur/' | relative_url }}">Continue to film and visual storytelling <span aria-hidden="true">↗</span></a>
  </section>
</div>

<div class="artist-gallery__lightbox" data-gallery-lightbox aria-hidden="true" role="dialog" aria-modal="true" aria-label="Photograph viewer">
  <button class="artist-gallery__lightbox-close" type="button" data-gallery-close aria-label="Close photograph viewer">Close <span aria-hidden="true">×</span></button>
  <button class="artist-gallery__lightbox-nav artist-gallery__lightbox-nav--prev" type="button" data-gallery-prev aria-label="Previous photograph">←</button>
  <div class="artist-gallery__lightbox-stage">
    <img data-gallery-lightbox-image alt="">
    <div class="artist-gallery__lightbox-copy">
      <p data-gallery-lightbox-meta></p>
      <h2 data-gallery-lightbox-caption></h2>
    </div>
  </div>
  <button class="artist-gallery__lightbox-nav artist-gallery__lightbox-nav--next" type="button" data-gallery-next aria-label="Next photograph">→</button>
</div>

<script src="{{ '/assets/js/gallery-artist.js' | relative_url }}" defer></script>
