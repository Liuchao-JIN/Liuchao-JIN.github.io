# Liuchao-JIN.github.io 全站内容审查与优化方案

> 审查对象：`D:\Documents\GitHub\Liuchao-JIN.github.io`  
> 对应站点：<https://liuchao-jin.github.io/>  
> 审查日期：2026-07-20  
> 文档性质：内容、信息架构、学术可信度、公开资料、SEO 与持续维护的综合方案  
> 说明：本文档只提出方案，没有直接修改网站页面；它应取代 2026-07-17 的通用版 `OPTIMIZATION_PLAN.md`，作为下一轮重构的主清单。

> **发布注意：**本文档包含内部审查结论和风险文件名，不应作为网站内容公开。实施时请把本文件和旧 `OPTIMIZATION_PLAN.md` 加入 Jekyll `exclude`，或移到不参与 GitHub Pages 构建的私有文档位置。

---

## 1. 结论先行

这个网站的优点不是“内容多”，而是已经积累了很强的研究、论文、教学、报告、获奖和个人经历素材。当前最大的障碍是：**这些素材仍以个人文件档案柜的方式存在，没有被组织成一个可信、聚焦、可持续维护的研究者品牌。**

目前有五个核心问题：

1. **主线不清楚**：首页只用一段文字罗列大量关键词，没有解释“你解决什么问题、如何解决、产生了什么结果、未来要去哪里”。
2. **多个数据源互相冲突**：论文页、`_publications` 集合、网页 CV、PDF CV、Talks 和 Teaching 页面更新时间不同，访客可能看到不同数量、不同年份和不同状态。
3. **公开范围失控**：课程资料、教材、试题答案、论文全文、成绩单、证书、会员证明和含个人信息的 PDF 被直接公开并进入 sitemap。
4. **大量过期内容仍以 Current / Upcoming 呈现**：这会直接损害学术主页最重要的品质——可信度。
5. **两个以上主页域名并存**：`liuchao-jin.github.io`、`jin-liuchao.github.io`、`jinliuchao.github.io` 同时出现在网站和附件中，分散搜索权重，也增加维护成本。

建议把网站重构为一个围绕以下叙事展开的专业主页：

> **Liuchao Jin develops AI-powered methods to design, manufacture, and validate intelligent matter—from architected materials and multimaterial 3D/4D printing to bio-inspired soft robotic systems.**

推荐的三大内容支柱为：

1. **AI-powered (meta)material discovery**
2. **Intelligent manufacturing systems**
3. **Bio-inspired soft robotic systems**

整个网站的每一页都应回答一个明确问题：

- Home：你是谁，为什么值得继续了解？
- About：你的研究逻辑和成长路径是什么？
- Research：你解决哪些问题，用什么方法，有什么证据？
- Publications：哪些成果最能代表你的贡献，完整记录在哪里？
- News：你最近在做什么？
- CV：你的正式经历和记录是什么？
- Contact：如何与你开展合作？

---

## 2. 本次审查覆盖范围

### 2.1 仓库规模

| 项目 | 当前数量 / 规模 | 审查含义 |
|---|---:|---|
| Git 跟踪文件 | 1,283 | 已超过普通静态个人主页的合理复杂度 |
| `_pages` | 24 | 其中包含正式页面、中文指南、模板演示页和重复页面 |
| `_posts` | 4 | 旅行、香港生活、申请指南、个人经历文章 |
| `_publications` | 13 | 仅 12 条能正常形成完整页面，且远少于主论文页 |
| `_talks` | 4 | 全部是 AcademicPages 模板假数据 |
| `_teaching` | 16 | 少于 Teaching 页面列出的 20 门课 |
| `_portfolio` | 3 | 主要是第三方课程资料和文献 PDF，不是专业研究作品集 |
| `files/` | 929 个文件 | 约 1.63 GB / 1.52 GiB |
| PDF | 437 | 当前全部进入生产 sitemap |
| 线上 sitemap | 501 个 URL | 其中 437 个为 PDF，占约 87.2% |

### 2.2 主要内容量

| 页面 | 当前内容特征 |
|---|---|
| 首页 `_pages/about.md` | 约 114 个空格分词词项，1 段介绍 + 1 张图片 |
| 网页 CV `_pages/cv.md` | 约 1,179 个词项、106 个链接，信息密度过高 |
| 主论文页 `_pages/mypublication.html` | 32 条手写成果、182 个链接、32 张图 |
| 第二论文页 `_pages/publications.md` | 仅渲染 13 条旧 collection 数据 |
| Talks | 23 个活动，但没有稳定区分“演讲、展示、参会” |
| Teaching | 20 门课程，多数仅链接 syllabus |
| Research Literature Review | 约 15,398 个词项、230 个链接、75 张图 |
| PhD Course Materials | 157 个链接 |
| Undergraduate Course Materials | 61 个链接 |

### 2.3 自动化内容检查结果

- 实际内容中约有 **588 个 HTTP 引用**，应改为 HTTPS 或站内相对路径。
- 已确认 **13 个站内失效引用**，另有 **21 个明确 404 的外部目标**。
- 实际公开内容约有 **144 张图片，其中 143 张没有有效 alt**。
- 线上 62 个主题页面中，标准 `meta description` 为 **0/62**，`og:image` 为 **0/62**。
- 60 个会生成页面的内容文件中，45 个没有独立 description/excerpt。
- 12 个路由出现多个 H1；147 个链接使用 `Web`、`PDF`、`Download here` 等泛化文字。
- 117 个相关 UTF-8 文本文件未发现明显乱码，这是当前内容层面的正面项。

---

## 3. P0：发布任何新版之前必须先处理

这一节不是“锦上添花”，而是上线前的最低风险门槛。

### 3.1 立即下线敏感证明材料

`_pages/cv.md` 直接链接了学位证、完整成绩单、排名证明、会员证明和 PDF CV。抽样检查确认，公开文件中出现了或可能出现：

- 学号、出生日期、性别；
- 完整成绩和课程记录；
- 签名、盖章、证书编号；
- 个人电话号码；
- 会员编号和有效期；
- 内部表格、会议链接、学生标识符。

建议：

1. 从当前部署和公共仓库中移除成绩单、排名证明、学位证、会员证书、身份证明、临时身份证、内部表格等文件。
2. 网页只保留文字事实，例如 `GPA 4.0/4.0 · WAM 96.29/100 · Ranked 1/79`，不直接链接扫描件。
3. 如果确实需要第三方核验，优先链接学校或颁奖机构的公开页面；没有公开核验页时，宁可不放证明附件。
4. 删除文件后要记住 Git 历史仍会保留旧对象；在做好备份、确认远端和协作者状态后，单独安排历史清理。
5. 需要时向搜索引擎申请缓存移除；`robots.txt` 和 `noindex` 不能代替访问控制。

### 3.2 立即暂停公开课程资料和第三方全文

当前公开内容包括：

- 商业教材和 solution manuals；
- 教师课件、作业、期中/期末考试及答案；
- candidacy examination 材料；
- 大量出版商论文 PDF、复制的摘要和图；
- 商业书籍全文；
- 由第三方电影片段组成的衍生视频描述。

重点位置：

- `_portfolio/PhD_Course_Materials.md`
- `_portfolio/Undergraduate_Course_Materials.md`
- `_portfolio/Research_Literature_Review.md`
- `_pages/amateur.md`
- `files/course_materials/`
- `files/essay/`
- `files/amateur/`

这些页面虽然没有出现在主导航，但 `_portfolio` 配置为 `output: true`，仍可被访问和索引。建议先全部停止公开，再逐项做权利审查。以后只恢复：

- 自己创作且愿意公开的内容；
- 已取得明确公开许可的内容；
- 有清晰开放许可证的内容；
- 论文的 DOI / publisher / 合法 open-access manuscript 链接。

仓库现有 MIT `LICENSE` 仅适用于相应软件代码，并不自动授权第三方教材、论文、试题和图片的再分发。

### 3.3 修复会直接损害学术可信度的论文错误

发布新版前至少修复：

| 位置 | 当前问题 | 建议 |
|---|---|---|
| `_pages/mypublication.html:132–137` 与 `260–265` | 两篇不同论文使用同一 DOI | 后一篇核对为 `10.1088/1361-665X/adcd19`；逐条对 DOI landing page |
| `_publications/jin2024machine.md` | article number 写成 `113086` | 应核对并改为 `102373` |
| `_publications/zhang2024high2.md` | YAML 因 `Poisson's` 单引号冲突而解析失败；PDF 指向另一篇论文 | 改用双引号或 block scalar，并修正 PDF |
| `_pages/mypublication.html` | “Machine learning powered inverse design…”作者顺序错误 | 以正式论文为准重新导入 |
| `zhai2023isogeometric` | 两个页面的期刊名称不一致 | 统一正式 bibliographic container |
| `jin2024spider` | 会议名称与出版容器混用 | 出版容器写 *MATEC Web of Conferences*，活动另列 |
| 多个 BibTeX | 作者拼写、姓/名顺序、加粗命令和正式标题不一致 | 从 DOI/Crossref/出版社导出，不再手写 |

主论文页目前有 32 条，`_publications` 有 13 条，PDF CV 有 34 条；审计后至少能识别 37 项 published/early-view work。**在完成逐条核验前，不建议在首页写死“37 publications”。**

### 3.4 删除 Under Review 目标期刊列表

当前论文页只列期刊名，不列标题，并公开 15 个 under-review 目标期刊。这样做有四个问题：

- 访客无法判断每个状态对应哪篇工作；
- 状态变化后极易过期；
- 可能与匿名评审或投稿沟通不协调；
- 目标期刊名本身不能证明成果。

建议规则：

- `Published / Early View`：完整引用 + DOI；
- `Accepted`：完整题目、作者、期刊、accepted 状态和日期；
- `Preprint`：有公开预印本链接时可列；
- `Under Review`：默认不列。确有公开预印本时，只写 `Preprint`，不写目标期刊。

### 3.5 清除模板假内容和公开开发页面

以下内容不应出现在正式主页或 sitemap：

- `_talks/2012...2014` 四个 `Relevant Topic in Your Field` 假条目；
- `_data/authors.yml` 中的 `Name Name` 和 `Name2 Name2`；
- `_pages/archive-layout-with-content.md`；
- `_pages/markdown.md`；
- `_pages/non-menu-page.md`；
- `_pages/talkmap.html` 及其中的假地点；
- `talkmap/map.html` 的 `Leaflet debug page`；
- `_drafts/post-draft.md` 的 lorem ipsum；
- 模板评论和空白测试页；
- 已发布的 `OPTIMIZATION_PLAN` 页面。

### 3.6 核实并统一当前身份

现在同时存在以下表述：

- 首页和 author profile：`PhD Candidate @ CUHK`；
- HTML CV：`2022–2026 (expected)`；
- 首页：仍“currently”在 SUSTech 与 Shenzhen University 访问；
- HTML CV：两段访问经历均在 2025-09 结束；
- PDF CV：又包含更多未同步经历。

在修改任何个人定位前，先做一个明确选择：

| 实际状态 | 推荐用语 |
|---|---|
| 学位尚未正式授予 | `PhD candidate in Mechanical and Automation Engineering at CUHK` 或 `Doctoral researcher at CUHK` |
| 已通过答辩但学位待授予 | `PhD candidate (degree expected [month year])`，不要提前写 `PhD` |
| 学位已经授予 | `Researcher in Intelligent Matter & Advanced Manufacturing`，About/CV 中写正式授予日期 |
| Fellowship 只是获奖、任职尚未开始 | 放入 Recognition，不写成 current postdoctoral position |
| 任职已正式开始 | 以正式机构、职位和起始月为准更新全站 |

用户此前已明确要求移除所有 Caltech 内容，因此本轮重构应从 HTML、PDF CV、博客、元数据和附件描述中一并删除 Caltech；不要因为旧 PDF CV 中存在该经历而重新引入。

### 3.7 更新高时效指南，或直接归档

中文指南包含申请费、入境申报、银行卡、电话卡、交通补贴、医疗和校园手续等高时效信息。至少已有明确过期例子：

- CUHK 申请指南仍写 HK$300；当前官方页面显示每个 programme 的申请费为 HK$500：[CUHK Graduate School Fees](https://www.gs.cuhk.edu.hk/admissions/scholarships-fees/fees)。
- 文章仍要求日常出入境健康申报；常规健康申报已自 2023-11-01 起取消：[中国政府网说明](https://english.www.gov.cn/news/202310/31/content_WS6540a78ec6d0868f4e8e0d15.html)。
- 交通补贴门槛仍写 HK$200；官方文件显示门槛自 2025-06-01 起调整为 HK$500：[公共交通费用补贴计划官方通知](https://www.ptfss.gov.hk/doc/change_threshold20250601.pdf)。

建议二选一：

1. **推荐**：把这些页面归档为 `Archived guide — not maintained`，只保留官方入口；
2. 如果继续维护：每页必须有 `Last verified`、`Next review`、官方来源、免责声明和维护责任人，且至少每半年复核一次。

不要再推荐把身份证明或申请材料上传至第三方在线 PDF 工具；优先建议本地离线处理。

---

## 4. 先决定唯一主站

当前至少出现三个域名形态：

- `liuchao-jin.github.io`：本仓库配置的域名；
- `jin-liuchao.github.io`：新版主页及 PDF CV 中出现；
- `jinliuchao.github.io`：旧课程链接中出现。

### 推荐方案

选择一个唯一的 professional homepage，并把其他域名做永久重定向或明确标记为 archive。不要维持三个内容近似、数据不同的个人主页。

| 选择 | 适用情况 | 行动 |
|---|---|---|
| 以新版 `jin-liuchao.github.io` 为主站 | 新版视觉、Research、News、Contact 和交互已经更成熟 | 本仓库只保留重定向或内容归档；所有 CV/Scholar/ORCID/GitHub 统一指向新版 |
| 以本仓库 `liuchao-jin.github.io` 为主站 | 希望保留 Jekyll/AcademicPages 和当前域名积累 | 把新版成熟内容迁回本仓库，另一个域名重定向到这里 |
| 两站分工 | 只有在确有长期维护能力时使用 | 一个明确为 professional homepage，另一个明确为 Chinese notes/archive；必须用不同标题、导航和 canonical |

本文后续方案假设：**本仓库仍要作为完整 professional homepage 重构**。如果最终选择新版为主站，本方案仍可作为内容清理和迁移清单使用。

---

## 5. 目标受众与内容任务

| 主要访客 | 他们最关心什么 | 网站应在何处回答 |
|---|---|---|
| 潜在导师、PI、招聘委员会 | 研究定位、独立贡献、未来方向、代表成果 | Home、Research、Selected Publications、CV |
| 潜在合作者 | 方法能力、设备/制造能力、可合作问题、联系方式 | Research、Contact |
| 学术同行和审稿人 | 完整且准确的论文、DOI、作者角色、代码/数据 | Publications |
| 会议组织者和媒体 | 简短 bio、头像、talk topics、荣誉、联系方式 | About、Talks、Contact |
| 学生和年轻研究者 | 研究解释、教学经验、合法公开资源 | Research、Teaching、Writing |
| 普通访客 | 你做的研究为什么重要 | Home、Research 的可视化和案例 |

内容设计原则：

1. 首先展示“贡献”，不是堆关键词。
2. 首先展示“证据”，不是堆形容词。
3. 完整信息放内页，首页只保留最高价值内容。
4. 每个事实只有一个权威数据源。
5. 每个时间敏感内容都有复核日期和过期机制。

---

## 6. 推荐信息架构

### 6.1 主导航

```text
Home · About · Research · Publications · News · CV · Contact
```

### 6.2 次级入口

以下内容不应抢占主导航：

- Teaching & Talks
- Writing / Notes
- Beyond Research
- Chinese Guides / 中文指南
- Privacy
- Sitemap

可放在页尾或 `More` 菜单中。

### 6.3 推荐 URL

| 内容 | 推荐 URL | 旧 URL 处理 |
|---|---|---|
| Home | `/` | 保留 |
| About | `/about/` | 当前首页 redirect 结构需拆分 |
| Research | `/research/` | 新建核心页 |
| Publications | `/publications/` | `/my_publication/` 301/redirect 到此 |
| News | `/news/` | 新建 |
| CV | `/cv/` | 保留，但重写 |
| Contact | `/contact/` | 新建 |
| Teaching | `/teaching/` | 保留 |
| Talks | `/talks/` | 保留 |
| Writing | `/writing/` | 取代 `Blog Posts` |
| Beyond Research | `/beyond-research/` | 取代 `/amateur/` |
| Chinese Guides | `/zh/guides/` | 合并重复指南；设置 `lang: zh-CN` |

---

## 7. Home：从“欢迎页”改成研究者入口

### 7.1 当前问题

`_pages/about.md` 只有一段介绍，存在：

- `Welcome to Liuchao Jin's Homepage!` 过于模板化；
- `Mechanical & Automation Engineering` 不如正式名称统一；
- `act as visiting scholars` 有语法和时态问题；
- 一次列出 3D/4D printing、smart materials、metamaterials、soft robotics、energy absorption、energy harvesting，缺少层级；
- 没有代表成果、最新动态、清晰 CTA 或研究方法；
- `Cover Image` 不能解释图片内容、来源和关联论文。

### 7.2 推荐首屏文案

页面标签：

```text
Liuchao Jin | Homepage
```

Eyebrow：

```text
AI · MECHANICS · MULTIMATERIAL MANUFACTURING
```

姓名：

```text
Liuchao Jin
```

身份行：

```text
Researcher in Intelligent Matter & Advanced Manufacturing
```

主标题：

```text
From AI-driven intelligent matter to autonomous scientific discovery
```

简介：

```text
I combine machine learning, computational mechanics, and multimaterial 3D/4D printing to design adaptive structures and soft robots. My long-term goal is to connect scientific reasoning, simulation, fabrication, and experiments in an evidence-grounded discovery loop.
```

CTA：

```text
Explore Research · View Publications · Download CV
```

可信度信息条：

```text
Doctoral research
CUHK · Hong Kong PhD Fellowship Scheme

Research experience
CUHK · SUSTech · Westlake · McGill · SCU

Recent recognition
2026 Best Presentation & Best Poster Awards
```

说明：用户此前明确要求首页 Research experience 不展示 Shenzhen University，因此首页应使用以上五项；完整 CV 是否保留 Shenzhen University 的历史经历，应在核实事实和展示意图后单独决定。

### 7.3 首页推荐顺序

1. Full-height Hero：一句话定位 + 主视觉 + 3 个 CTA。
2. Three Research Pillars：三张卡片，每张一句问题、一句方法、一项代表证据。
3. Selected Research：4–6 项最能体现个人贡献的成果。
4. Research Loop：Reason → Design → Simulate → Fabricate → Validate → Learn。
5. Selected News：最新 4–6 条。
6. Recognition：只放最近且重要的 4–6 项。
7. Contact CTA：合作、邀请报告、研究访问。

不要在首页加入完整论文列表、完整奖项列表、所有社交媒体或所有访问经历。

---

## 8. About：建立完整、可信且有高度的个人叙事

### 8.1 页面目标

About 不是 CV 的散文化复制，而是解释：

- 你的核心科学问题是什么；
- 你的方法为什么不同；
- 你的研究如何从材料走向制造和机器人；
- 你未来想构建什么系统；
- 哪些经历支撑这个方向。

### 8.2 推荐长版 Profile 草稿

以下英文草稿可作为 About 的主体。上线前必须把方括号内状态核实并替换：

```text
Liuchao Jin is a researcher at the intersection of artificial intelligence, computational mechanics, advanced manufacturing, and embodied intelligence. [If the degree has not yet been conferred: He is completing his doctoral research in Mechanical and Automation Engineering at The Chinese University of Hong Kong under Professor Wei-Hsin Liao as a Hong Kong PhD Fellowship awardee.] His academic foundation was established at the Sichuan University–Pittsburgh Institute, where he completed a BEng in Mechanical Engineering with a GPA of 4.0/4.0, a weighted average mark of 96.29/100, and a rank of 1/79.

His research begins with a central question: given a desired physical behavior, can an algorithm discover a material architecture that is not only mechanically effective but also manufacturable and experimentally verifiable? To answer this question, he integrates machine-learning surrogate models, mechanics-informed optimization, finite-element simulation, multimaterial 3D/4D printing, and physical experiments. This workflow transforms target shapes, strain fields, stiffness profiles, and energy-absorption requirements into realizable material systems.

Across architected materials, programmable structures, intelligent additive manufacturing, and bio-inspired robotics, his work treats algorithms, geometry, materials, fabrication processes, and validation as one coupled system. Rather than optimizing a digital design in isolation, he focuses on closing the loop between computational intent and physical performance.

His longer-term vision is to develop evidence-grounded systems for autonomous scientific discovery—systems that can coordinate reasoning, hypothesis generation, simulation, fabrication, experiments, and knowledge refinement while preserving traceability, uncertainty awareness, and meaningful human oversight.
```

### 8.3 About 页面建议模块

1. **Professional Profile**：上面的 4 段正文。
2. **Central Question**：一句大标题 + 一段解释。
3. **Research Principles**：
   - Physics-grounded intelligence
   - Manufacturability from the start
   - Evidence through validation
4. **Education & Research Journey**：时间线，不用证书扫描件。
5. **Selected Recognition**：重点荣誉。
6. **Downloadable Bio**：50-word / 100-word / full bio，方便会议组织者使用。
7. **Portrait & media kit**：一张正式头像，给出摄影来源和允许用途。

### 8.4 Recognition 推荐内容

优先展示：

- Presidential Global Impact Postdoctoral Fellowship Scheme 2025/26
- Best Presentation Award, 5th 4D Materials Design and Additive Manufacturing Conference, 2026
- Best Poster Award, ICAST 2026
- Reaching Out Award, HKSAR Government Scholarship Fund, 2026
- PhD International Mobility for Partnerships and Collaborations Award, 2025
- Outstanding Students Award, The Chinese University of Hong Kong, 2025
- Best Paper Award, IEEE ICUS, 2022
- Hong Kong PhD Fellowship Scheme
- National Scholarships, 2019–2021（正式年份需按证书复核）

奖项统一格式：

```text
Official Award Name — Issuing Body, Month Year
One sentence explaining selectivity or purpose, only when a source supports it.
```

不要在 Recognition 中放奖金额、证书编号或扫描件。Fellowship 如果尚未开始任职，应继续作为 award，而不是 current position。

---

## 9. Research：把研究方向变成可理解的科学故事

### 9.1 页面总标题

```text
From intelligent matter to autonomous scientific discovery
```

简介：

```text
My current work co-designs algorithms, material architectures, and manufacturing processes for adaptive structures and soft robots. My long-term vision is to connect scientific reasoning, simulation, fabrication, and experiments in a governed, evidence-grounded loop.
```

### 9.2 三大研究方向标准文案

#### 01 · AI-powered (meta)material discovery

```text
Learning structure–response relationships and using physics-informed optimization to discover manufacturable architectures that realize target shapes, strain fields, stiffness profiles, or energy-absorption behavior.
```

Tags：

```text
Inverse design · Machine learning · Multiscale mechanics
```

#### 02 · Intelligent multimaterial manufacturing

```text
Developing digitally controlled 3D/4D printing processes for multimaterial structures, including metal–polymer systems, ceramics, lightweight tooling, and programmable digital materials.
```

Tags：

```text
3D/4D printing · Process intelligence · Digital materials
```

#### 03 · Bio-inspired soft robotic systems

```text
Combining origami-inspired geometry, hard–soft material coupling, and adaptive actuation to create compliant mechanisms that transform material intelligence into embodied action.
```

Tags：

```text
Soft robotics · Origami · Embodied intelligence
```

### 9.3 代表研究不要只放摘要

每项案例统一使用：

```text
Challenge — What physical behavior or manufacturing constraint made the problem difficult?
Approach — What did you personally develop or integrate?
Outcome — What was demonstrated, measured, fabricated, or validated?
Evidence — DOI · Code · Data · Video, when available and authorized
```

建议优先准备四个案例：

1. Inverse design of strain fields in hierarchical architectures.
2. Inverse design of stimulus-driven shape change.
3. Multimaterial printing of low-melting-point alloy/polymer systems.
4. Origami-inspired hard–soft robotic grippers.

每个案例要写清个人贡献，不要直接复制论文 abstract。

### 9.4 长期愿景

把“AI for Science”从口号变成可验证的 research programme：

```text
Reason → Hypothesize → Design → Simulate → Fabricate → Challenge → Learn
```

并解释三个治理原则：

- Traceable evidence
- Explicit uncertainty
- Human oversight at consequential decisions

### 9.5 互动内容

可在 Research Pillars 后加入 `Living Research Matter` 入口：

```text
Experience the research system
Watch particles evolve from material discovery to multimaterial fabrication and embodied intelligence.
```

互动页必须仍有普通 HTML 文案和链接作为降级内容；不能让 Canvas 成为唯一的信息载体。

---

## 10. Publications：从三套冲突数据改成一个可信系统

### 10.1 当前冲突

| 数据源 | 数量 / 状态 | 问题 |
|---|---:|---|
| `_pages/mypublication.html` | 32 条 | 主导航使用；全手写；有 DOI/作者/年份错误 |
| `_publications/` | 13 条，约 12 条完整 | CV 和第二论文页使用；缺少近两年大量成果 |
| PDF CV | 34 条 | 与网页不同；仍有 under-review 列表和标题错误 |
| 文件夹 PDF/BibTeX | 更多 | 文件存在不等于网页记录完整或正确 |

### 10.2 唯一数据源

建议使用 `_data/publications.yml`、BibTeX/JSON，或彻底重建 `_publications`。每篇记录至少包含：

```yaml
slug:
title:
authors:
  - name:
    self: true|false
    equal_contribution: true|false
    corresponding: true|false
status: published|early-view|accepted|preprint
type: journal-article|conference-paper|review
online_date:
issue_date:
venue:
volume:
issue:
pages_or_article_number:
doi:
publisher_url:
pdf:
pdf_rights: author-manuscript|open-access|not-hosted
bibtex:
thumbnail:
thumbnail_alt:
selected: true|false
contribution:
awards:
```

由同一数据源生成：

- 首页 Selected Publications；
- `/publications/` 完整列表；
- CV 的 Selected Publications；
- 每篇独立页面；
- BibTeX 下载和 JSON-LD。

### 10.3 年份规则

很多论文存在 online-first 年份和 final issue 年份冲突。建议同时保存：

- `online_date`：首次正式在线日期；
- `issue_date`：卷期正式日期；
- 页面分组：以 final issue year 为主；
- early view 尚无卷期时：暂按 online year，并标记 `Early View`。

### 10.4 页面结构

1. 顶部一句话说明研究覆盖范围。
2. 6–8 项 Selected Publications，每项有贡献说明。
3. 筛选器：Year / Research Pillar / Type / Role。
4. 完整列表，按年份倒序。
5. Scholar / ORCID 入口。

推荐的精选成果类型：

- 2026 functional-material inverse-design review；
- 2025 strain-field inverse design；
- 2025 soft-robot FEA/ML/digital-twin review；
- 2024 arbitrary-shape 4D-printing inverse design；
- 2024 big-data/ML/digital-twin additive-manufacturing review；
- 2025 INPR connector work；
- low-melting-point alloy additive manufacturing；
- 2022 UAV controller with Best Paper Award。

最终选择应按“个人贡献 + 三大方向覆盖 + 代表性证据”决定，而不是只看期刊名。

### 10.5 链接和标记规范

- 使用 `https://doi.org/...`，标签写 `DOI`，不用模糊的 `Web`。
- PDF 只在有权公开时托管；否则链接 DOI 或合法 repository。
- `*`、`†` 只在页面顶部定义一次。
- ESI Highly Cited、Cover Paper、Top Downloaded 等作为有日期和来源的 badge，不放进论文标题。
- BibTeX 从 DOI/出版社生成；不要在 BibTeX 内用 `\textbf{}` 标记本人。

### 10.6 2026 状态复核

实施时优先检查当前 Accepted 区的三项是否已经转为 Early View/Published，并通过 DOI 页面补齐题目、作者和日期。审计发现的候选 DOI 包括：

- `10.1002/adma.74024`
- `10.1016/j.ijmecsci.2026.111900`
- `10.1016/j.tws.2026.115252`

这些信息具有时间敏感性，上线前必须再以出版社页面核验。

---

## 11. News：让网站体现活跃度，但不制造状态债务

### 11.1 推荐字段

```yaml
date: 2026-07-10
type: award|paper|talk|visit|service
title:
summary:
link:
image:
image_alt:
expires_at:
featured: true|false
```

### 11.2 推荐首批内容

在核实正式名称和日期后，可包括：

- Best Presentation Award at 4DMDA 2026, Dublin.
- Best Poster Award at ICAST 2026.
- Presidential Global Impact Postdoctoral Fellowship Scheme 2025/26.
- Reaching Out Award, HKSAR Government Scholarship Fund.
- 2026 papers moved to Early View / Published.
- Selected invited talk or conference presentation.

### 11.3 内容规范

- 最新在上；首页展示 4–6 条，完整页显示全部。
- 每条必须说明“发生了什么”，不要只写期刊或活动名称。
- 不使用 `Upcoming` 的永久手写分类；由日期自动计算。
- 过期活动自动移动到 Past。
- 每条尽量链接官方会议、论文 DOI 或机构新闻，不链接敏感证书。

---

## 12. CV：保留完整性，但突出决策信息

### 12.1 推荐顺序

1. Current position and contact
2. Research interests
3. Education
4. Research and professional appointments
5. Selected publications
6. Patents
7. Selected honors
8. Teaching
9. Selected professional service

### 12.2 应删除或压缩

- 删除成绩单、排名、学位证、会员证书的直接链接。
- 4 页 CV 只保留 10–15 项 selected publications；完整列表链接 Publications/Scholar。
- 删除 under-review 目标期刊和大段未公开 manuscript 列表。
- 约 29 家 reviewer 期刊压缩为 `Selected peer-review service`；如有 ORCID/Web of Science 记录可链接。
- 奖项只留 5–8 个最重要项目；完整清单可以放 About/Recognition 内页。
- 奖金额度、证书编号、献血等与岗位评估关联弱的信息移出 professional CV。

### 12.3 必须复核的事实

- PhD 是否已经正式授予；
- SUSTech、Shenzhen、Westlake、McGill/Mitacs 的正式起止时间和身份；
- Teaching Assistant 是否在 2025-07 结束；
- HKYSA、HKQTMA 等会员有效期是否已结束；
- PDF CV 和 HTML CV 中的论文题目、作者顺序与年份；
- 所有 Caltech 相关表述按用户要求删除；
- Presidential Global Impact Postdoctoral Fellowship Scheme 2025/26 应放 Recognition，除非正式任职已开始。

---

## 13. Talks：把参会清单改成学术传播记录

### 13.1 当前问题

- 2026-03 的会议在 2026-07 仍列为 Upcoming；
- 多数项目只写会议和城市，没有 presentation title 或 role；
- 宁波与 Plymouth、厦门与成都出现日期重叠，需要说明线上、部分参会或修正日期；
- `Shangdong` 拼写错误；
- 至少一个 talk PDF 缺失；
- `_talks` 内仍有四条模板假数据。

### 13.2 推荐分类

1. Invited Talks
2. Conference Presentations — Oral
3. Conference Presentations — Poster
4. Panels / Academic Forums
5. Conference Attendance / Service（可选，默认不放主页面）

### 13.3 每条模板

```text
Presentation title
Role · Event · Location / Online · Date
One sentence on the topic or outcome
Abstract · Slides · Recording, when public and authorized
```

示例：

```text
Machine-learning-enabled design of hierarchical architectures for 4D printing
Invited talk · GAMES Webinar 383 · Online · 30 Oct 2025
Recording
```

参会但没有演讲的活动不要放在 `Talks & Presentations` 中冒充 presentation。

---

## 14. Teaching：从 syllabus 清单改成教学能力证据

### 14.1 当前问题

- Spring 2025 仍标为 Current；
- Teaching 页面有 20 门，collection 只有 16 门；
- MAEG5735 的 linked PDF 标题与页面课程名冲突；
- 两门 4xxx 课程被错误写成 graduate course；
- 多个 syllabus 没有体现本人角色；
- 部分 PDF 暴露学生标识、内部群组或联系信息。

### 14.2 推荐表达

```text
Teaching Assistant, The Chinese University of Hong Kong — Aug 2022–Jul 2025

MAEG5735 Applied Computational Intelligence — Spring 2025
Role: [tutorials / labs / grading / office hours — verify]
Topics supported: [verify]
Teaching contribution or outcome: [verify]
Official course catalog · Original teaching material, if authorized
```

重复开设可合并：

```text
MAEG5725 Control and Industrial Automation — Spring & Fall 2024
```

### 14.3 教学资料原则

- 链接官方 course catalog；
- 只公开本人原创、已脱敏并获许可的教学材料；
- 不公开 syllabus 作为“任教证明”，除非它确实标明本人角色且允许公开；
- 不公开学生名单、学号、群组、会议链接、试题和答案；
- 自己的讲义应注明许可证，例如 `CC BY 4.0`。

---

## 15. Writing 与中文指南

### 15.1 重新分区

建议把 `Blog Posts` 改名为 `Writing`，分为：

- Research perspectives
- Personal essays
- Field notes
- Chinese guides / 中文指南（独立语言区）

### 15.2 香港生活与申请指南

`_posts/2023-08-08-hkinfo.md` 与多个 `_pages` 指南重复，并包含高时效流程。推荐合并为一个 canonical 页面，页面顶部固定显示：

```text
Archived / Last verified: 20 Jul 2026
This guide reflects the author’s experience at the time. Policies, fees, and procedures change; official CUHK and HKSAR sources take precedence.
```

每条程序性信息必须链接官方来源。银行、返现、电话套餐、价格和医疗体验属于个人经验，应明确写 `personal experience, not an endorsement`，或直接删除。

`“这边随便填”` 应改为：

```text
请根据实际情况如实填写；申请人应对提交信息的准确性负责。
```

### 15.3 旅行内容

`_pages/travel.md` 与 2023 博文相似度约 95.7%，只保留一份。推荐开头：

```text
In June 2023, I visited Jinan and Qingdao as part of the Hong Kong Young Scholars Shandong tour. This photo essay records the research exchanges and cultural sites that shaped the visit, including [two verified examples].
```

需要：

- 修正 `/files/trip/` 与实际 `/files/Trip/` 的大小写冲突；
- 每张照片写 alt 和 caption；
- 删除 EXIF 中的 GPS、拍摄时间和设备信息；
- 确认照片中可识别人物的公开授权。

### 15.4 个人经历文章

`From Jin Village to Global Stage` 的优点是有人格和成长背景，但应优化：

- 将祖源表述明确为 `family oral history`，除非有可靠史料；
- 删除关于“中国父母”的宽泛 stereotype；
- 合并两段重复的 survival / dignity / education 论述；
- 从“资源匮乏叙事”转向具体行动：如何影响你的研究、教学、导师关系或支持年轻学者；
- 按用户要求删除 Caltech；
- 修正 permalink 中 `/2025/30/` 把日期误作月份的问题。

可考虑新标题：

```text
From Jin Village to Research: Education, Mobility, and Scientific Belonging
```

---

## 16. Beyond Research：不要再叫 Amateur

`Amateur` 容易让高质量个人创作显得次要或不专业。建议改成：

```text
Beyond Research
```

推荐简介：

```text
Outside the lab, I read about gender, mythology, and speculative fiction, and I produce short films and campus news videos. These projects sharpen how I think about narrative, technology, and whose experiences are represented in scientific culture.
```

内容策略：

- 删除商业书籍 PDF，改为出版社、图书馆或书目链接；
- 精选 3–5 个最好的视频，而不是列出所有视频；
- 每项说明日期、本人角色、合作者、作品背景和版权/素材来源；
- remix 内容需要明确素材权利和创作边界；
- 社交平台只保留真正活跃且与此板块相关的入口。

---

## 17. Portfolio / Resources：先下线，再重建为原创资源

当前 Portfolio 更像未经治理的文件目录，不适合作为职业作品集。推荐重建为：

### 17.1 Research Notes

把 `Research Literature Review` 从复制摘要和托管 PDF 改成原创 annotated research map：

- 研究问题；
- 领域脉络；
- 每篇论文 2–3 句原创评价；
- 方法对比表；
- 开放问题；
- DOI 链接；
- `Last updated` 和版本记录。

### 17.2 Teaching Resources

只保留本人原创且授权公开的：

- tutorial notes；
- example code；
- interactive demonstrations；
- reading lists；
- assignment templates without answers。

### 17.3 Research Artifacts

大型代码和数据不应继续塞入 GitHub Pages 仓库：

- 代码：独立 GitHub repository；
- 数据：Zenodo / Figshare / OSF 等带版本和 DOI 的平台；
- 视频：合适的视频托管；
- 站点只保留摘要、缩略图和稳定链接。

---

## 18. Contact：完整但克制

### 18.1 推荐内容

```text
Liuchao Jin
Researcher in Intelligent Matter & Advanced Manufacturing

Email
liuchao.jin@link.cuhk.edu.hk

Office
Room 201, William M.W. Mong Engineering Building (ERB)
The Chinese University of Hong Kong
Shatin, New Territories, Hong Kong SAR

Best for
Research collaboration, invited talks, academic visits, and scientific exchange.
```

### 18.2 公开位置的平衡

用户此前希望显示 ERB 201 和英文 Google Map。推荐只在 Contact 页面显示完整办公地址和英文地图，不要在所有页面的 author profile 重复精确房间号。上线前确认该办公室和机构身份仍有效。

### 18.3 个人资料链接

核心链接只保留：

- Google Scholar
- ORCID (`0000-0001-6204-1922`)
- ResearchGate（可选）
- LinkedIn
- GitHub
- Email
- Book with me（确认链接仍有效后）

Facebook、Instagram、X/Twitter、YouTube、Pinterest 如与专业主页关系弱，应移到 Beyond Research 或删除。Website 自链接没有价值。

---

## 19. Footer、Privacy、Sitemap 与 404

### 19.1 Footer

所有页面使用同一模块化 footer，内容建议为：

- 一句研究定位；
- About / Research / Publications / News / Contact；
- Email / Scholar / ORCID / LinkedIn / GitHub；
- CUHK 标识只在符合品牌规范并获准使用时展示；
- Copyright；
- Privacy；
- 可选的访问量。

### 19.2 访问量与 ClustrMaps

`author-profile.html` 仍嵌入 ClustrMaps，而此前服务已经不稳定。建议移除该脚本。如果保留历史访问量：

```text
Global visits 82,456+
Legacy count carried over from ClustrMaps; new counting started [date].
```

不要把历史基线和新统计悄悄相加。新统计工具应：

- 尽量少收集数据；
- 不展示精确访客位置；
- 在 Privacy 中如实说明；
- 加载失败时不影响页面内容。

### 19.3 Privacy

当前 2016 模板错误声称使用 Disqus、广告商和 Google Analytics，同时没有说明实际 ClustrMaps。应重写为基于实际服务的简短隐私说明：

- 托管方；
- 是否使用 analytics/counter；
- 是否设置 cookie；
- 第三方嵌入（Google Maps、YouTube 等）；
- 联系方式；
- 最后更新日期。

### 19.4 Sitemap

XML sitemap 只收录真正公开的内容页面，不收录：

- 证书、成绩单和内部文件；
- 课程资料；
- 测试页、模板页、优化方案；
- 直接 PDF 资产，除非有明确理由；
- 重复页面和归档页面。

HTML sitemap 应是人工策划的内容目录，不要遍历所有 collections 和 files。

### 19.5 404

删除旧 Google fixURL 脚本。提供：

```text
The page may have moved.
Home · Research · Publications · Contact
```

---

## 20. 统一文风与术语

### 20.1 人称

- Home、Research、Contact、Writing：第一人称，直接、有行动感。
- 可复用会议 bio：第三人称。
- CV：省略主语的正式条目。

### 20.2 推荐术语

| 不统一写法 | 推荐写法 |
|---|---|
| Ph.D. / Ph.D | PhD |
| BE | BEng |
| Mechanical & Automation Engineering | Mechanical and Automation Engineering |
| Sichuan University - Pittsburgh Institute | Sichuan University–Pittsburgh Institute |
| meta-materials / meta-structures | metamaterials / architected materials |
| machine learning enabled | machine-learning-enabled |
| multi-material | multimaterial（除非正式论文标题不同） |
| bioinspired | bio-inspired |
| hard-soft | hard–soft |
| 2022-2026 | 2022–2026 |

### 20.3 写作规则

- 每段只承担一个逻辑任务。
- 先写结果，再写方法背景。
- 少用 `novel`、`first`、`world-leading`；除非有正式来源支持。
- 不复制论文摘要；改写为访客能理解的 contribution statement。
- 论文标题严格使用出版社版本，包括英式/美式拼写。
- 日期统一为 `Jul 2026` 或 `20 Jul 2026`。
- 奖项写官方名称、颁发机构和年份。
- 所有 Current / Upcoming / Under Review 内容必须有自动状态或复核日期。

---

## 21. 图片、视频和附件内容标准

### 21.1 图片

每张内容图必须有：

- 描述性 alt；
- caption；
- 来源或作者；
- 对应论文/项目链接；
- 公开授权状态；
- 去除 EXIF；
- 合理尺寸和压缩格式。

示例：

```text
Alt: Predicted and experimentally measured strain fields for three target patterns in a hierarchical multimaterial architecture.

Caption: Machine-learning-powered inverse design connects target strain fields to manufacturable material allocations. Adapted from [paper title], 2025.
```

装饰图使用空 `alt=""`；不能用 `Cover Image`、`Image 1`。

### 21.2 论文图

- 核实 publisher reuse policy；
- 优先使用作者生成的重新绘制示意图；
- 不直接裁剪出版 PDF 作为缩略图；
- 每个 Selected Work 的图片保持相同画布比例。

### 21.3 视频与地图

- 视频提供 poster、字幕/文字说明、暂停控制和 reduced-motion 降级；
- 背景视频不能承载唯一信息；
- Google Map 使用英文界面，iframe 有 `title`；
- 第三方嵌入在 Privacy 中披露。

---

## 22. SEO、分享和结构化数据

### 22.1 全站基础元数据

建议：

```yaml
title: "Liuchao Jin | Homepage"
description: "Liuchao Jin is a researcher in AI-driven material design, computational mechanics, multimaterial additive manufacturing, and bio-inspired soft robotics."
locale: "en-US"
timezone: "Asia/Hong_Kong"
future: false
```

每页至少：

```yaml
title:
description:       # 120–160 English characters; not the full abstract
lang: en           # or zh-CN
permalink:
last_modified_at:
sitemap: true
image:
image_alt:
```

### 22.2 分享预览

- 一张统一的 1200×630 Open Graph 图；
- 每个 Research case / News item 可有独立图；
- `og:title`、`og:description`、`og:image`；
- Twitter Card；
- canonical 指向唯一主站。

### 22.3 Person JSON-LD

至少包含：

- name；
- url；
- image；
- affiliation（必须是当前事实）；
- alumniOf；
- sameAs：Scholar、ORCID、LinkedIn、GitHub；
- knowsAbout：三大研究方向。

不要把已经结束的 visiting appointment 写成 current affiliation。

### 22.4 现有配置需修复

- `title: "Jin's Homepage"` 过于泛化；
- `description: "personal description"` 是占位符；
- `future: true` 有误发未来内容风险；
- `timezone: America/Los_Angeles` 不符合香港内容；
- Universal Analytics 已过时且 tracking ID 为空；
- manifest 名称仍是 `Minimal Mistakes`；
- favicon 和 Android 图标有 7 个缺失资源；
- MathJax 2.7.4 全站加载，应只在需要公式的页面加载或升级。

---

## 23. 可访问性与内容语义

最低要求：

- 每页只有一个 H1；
- 中文页真正输出 `<html lang="zh-CN">`；
- 所有内容图片有 alt；
- 图表有 caption；
- 表格有 `<caption>` 和 `<th>`；
- 菜单按钮有可访问名称；
- iframe 有 title；
- 链接文字说明目标，例如 `Download the paper PDF`，不用 `Web` / `here`；
- 颜色不是唯一的信息表达；
- 键盘 focus 可见；
- 动画尊重 `prefers-reduced-motion`；
- 中文阅读时间不能继续按英文空格分词。

论文页目前用大量 table 和非法空 `<tr>` 做间距，应改为语义化 article/card 列表。

---

## 24. 内容数据架构

### 24.1 唯一数据源矩阵

| 内容 | 唯一数据源 | 生成位置 | 维护频率 |
|---|---|---|---|
| Biography / identity | `_data/profile.yml` | Home、About、Contact、JSON-LD | 身份变化时 |
| Research pillars | `_data/research.yml` | Home、Research、interactive page | 项目方向变化时 |
| Publications | `_data/publications.yml` 或规范化 `_publications` | Publications、CV、Home | 每篇状态变化时 |
| News | `_data/news.yml` | Home、News | 每月 |
| Awards | `_data/awards.yml` | About、CV、News | 获奖时 |
| Talks | `_talks` 规范记录 | Talks、CV、News | 活动后 7 天内 |
| Teaching | `_teaching` 规范记录 | Teaching、CV | 每学期 |
| Social/contact | `_data/profile.yml` | Header、Footer、Contact | 每季度 |

### 24.2 内容状态字段

对容易过期的项目增加：

```yaml
status:
start_date:
end_date:
last_verified:
review_by:
expires_at:
source_url:
```

页面上不再手写 `Current` 和 `Upcoming`，由日期计算。

### 24.3 CI 内容守门

每次 push 自动检查：

- YAML/front matter 严格解析；
- Linux 大小写路径；
- 内部链接和外链；
- 缺失 alt；
- 多 H1 和无标题 iframe；
- 禁止 `personal description`、`Name Name`、`Relevant Topic` 等占位词；
- 禁止未批准目录进入构建；
- 禁止新的 `http://Liuchao-JIN.github.io`；
- DOI 重复和必填字段；
- publication count 在各页面一致。

---

## 25. 分阶段实施路线

### P0：同日完成，先止损

- [ ] 备份仓库和远端状态。
- [ ] 暂停发布 `_portfolio` 和高风险附件。
- [ ] 下线成绩单、学位证、排名、会员证书、内部表格和含学生信息文件。
- [ ] 移除商业书籍、教材、试题答案、未授权论文全文。
- [ ] 删除 Under Review 目标期刊列表。
- [ ] 修复已知 DOI、YAML、PDF 指向和 publication 元数据错误。
- [ ] 删除模板 talks、虚构作者、开发/测试页面。
- [ ] 移除 Caltech 内容。
- [ ] 确认当前身份和唯一主域名。
- [ ] 移除 ClustrMaps，重写实际 Privacy。

### P1：1–3 天，重建内容骨架

- [ ] 新建 profile/research/publications/news/awards 唯一数据源。
- [ ] 重写 Home、About、Research、Contact。
- [ ] 合并两套 Publications 页面，设置 redirect。
- [ ] 重写 CV、Talks、Teaching。
- [ ] 新建 News。
- [ ] 重构主导航和模块化 footer。
- [ ] 把中文指南移入独立语言区或归档。

### P2：3–7 天，完成质量治理

- [ ] 逐条核验全部论文 DOI、作者、标题、年份、卷期、author role。
- [ ] 处理 588 个 HTTP 引用和已知失效链接。
- [ ] 修复路径大小写和旧域名。
- [ ] 为全部内容图补 alt、caption、credit，去除 EXIF。
- [ ] 更新 SEO、Open Graph、JSON-LD、favicon 和 manifest。
- [ ] 清理 sitemap、feed、404、重复页面。
- [ ] 建立 CI 内容检查。

### P3：持续优化

- [ ] 把 Research Literature Review 改成原创 annotated research map。
- [ ] 发布合法、原创的 teaching resources。
- [ ] 为 selected works 增加代码、数据、视频和 contribution statements。
- [ ] 将 Living Research Matter 作为 Research 的增强体验。
- [ ] 每季度做事实、链接、隐私、版权和状态审计。

---

## 26. 验收标准

新版上线前应满足：

### 内容与事实

- [ ] 全站只有一个 current role 和一套日期。
- [ ] PhD 状态已经核实。
- [ ] 所有面向访客的页面、CV 和公开附件中，Caltech 检索结果为 0（内部审计文档除外）。
- [ ] 首页 Research experience 为 `CUHK · SUSTech · Westlake · McGill · SCU`。
- [ ] 每项奖项使用正式名称和年份。
- [ ] Publications、CV、Home 的论文数量来自同一数据源。
- [ ] 不公开目标期刊形式的 Under Review 列表。

### 隐私与权利

- [ ] 生产站无法访问成绩单、证书、会员证明、身份证明和内部表格。
- [ ] 不再公开教材、solution manuals、考试答案和未授权论文 PDF。
- [ ] 公共照片无 GPS EXIF。
- [ ] 可识别人物和第三方素材有授权或已删除。
- [ ] Privacy 与实际 tracking/embed 完全一致。

### 信息架构

- [ ] 主导航最多 7 项。
- [ ] 每个主页面有明确任务和 CTA。
- [ ] 只有一个 Publications URL。
- [ ] 重复旅行/指南页面已合并或 redirect。
- [ ] 模板页、测试页和优化方案不进入 sitemap。

### 可访问性与 SEO

- [ ] 每页一个 H1。
- [ ] 所有内容图有有效 alt。
- [ ] 中文页输出 `lang="zh-CN"`。
- [ ] 无泛化 `Web` / `here` 链接。
- [ ] 每页有独立 description、canonical 和分享图。
- [ ] 无 HTTP mixed-content 资源。
- [ ] 内部链接在 Linux 大小写环境全部通过。

### 性能与维护

- [ ] 不再把整个 `files/` 目录当公共附件库。
- [ ] 大型数据/代码迁至适合的平台。
- [ ] CI 能阻止 YAML 错误、坏链接、占位内容和未批准附件。
- [ ] News、Talks、Teaching 和 publication status 有更新责任人。

---

## 27. 持续维护机制

### 每次新增论文

1. 从 DOI/出版社导入正式 metadata。
2. 填 contribution statement。
3. 检查 PDF 权利。
4. 准备一致比例缩略图和 alt。
5. 更新唯一数据源；让 Home/CV/Publication 自动同步。

### 每月

- 更新 News；
- 把过期 Upcoming 自动归档；
- 检查首页身份和 CTA；
- 跑内部链接和 YAML 检查。

### 每季度

- 外链检查；
- 论文状态、引用 badge 和奖项状态核验；
- 社交/预约/联系链接核验；
- 图片授权和公开附件审查；
- 搜索结果与 sitemap 抽查。

### 每年

- 重写 100-word bio 和研究定位；
- 重新选择 Selected Publications；
- 归档过时指南；
- 更新正式头像、OG 图和 CV；
- 审核 Privacy、analytics 和第三方嵌入。

---

## 28. 开始实施前需要本人确认的决策

1. 哪个域名是唯一 professional homepage？
2. PhD 是否已经正式授予？正式日期是什么？
3. 当前正式机构和职位是什么？
4. Presidential Global Impact Postdoctoral Fellowship 是否只是 award，还是任职已开始？
5. Shenzhen University 是否只从首页删除，还是也从完整 CV 删除？
6. 哪些论文 PDF 有明确公开权利？
7. 哪些教学材料是本人原创且允许公开？
8. ERB 201 是否仍是有效且愿意公开的办公地址？
9. Book with me 链接是否仍有效？
10. 中文指南是继续维护，还是统一归档？
11. 访问量 82,456 是只作为历史纪念显示，还是要接入新的统计？
12. Beyond Research 中哪些视频最能代表本人创作？

---

## 29. 推荐的最终内容交付物

完成本方案后，站点应形成以下可维护成果：

- 1 套唯一 profile 数据；
- 1 个聚焦的 Home；
- 1 个完整 About；
- 1 个三支柱 Research 页面；
- 1 套唯一 publication 数据和 1 个 canonical Publications URL；
- 1 个可自动归档的 News 系统；
- 1 份精简网页 CV + 1 份一致的 PDF CV；
- 规范的 Talks、Teaching、Contact；
- 合法、原创的 Resources；
- 清晰分区的 Writing / Chinese Guides / Beyond Research；
- 模块化 header/footer；
- 准确的 Privacy、Sitemap、404、SEO 和 JSON-LD；
- 一套自动内容质量检查。

最终目标不是把所有现有内容都“展示得更漂亮”，而是让每项公开内容都有明确价值、可信来源、正确状态、适当权限和稳定维护方式。这样访问者在几十秒内能理解你的研究主线，在几分钟内能验证你的贡献，并且不会因为旧页面、错误引用或敏感附件而降低对整个主页的信任。

---

## 30. 现有内容逐文件处置表

这一表格用于确保“所有内容都被审查”，实施时可直接逐行勾选。

### 30.1 `_pages`

| 文件 | 当前作用 | 建议处置 | 新内容方向 |
|---|---|---|---|
| `404.md` | 旧 404 + Google fixURL | 重写 | 简短说明 + Home/Research/Publications/Contact |
| `about.md` | 既是首页又充当 About | 拆分 | `/` 做 Home；`/about/` 做长版 Profile |
| `amateur.md` | 阅读和视频 | 改名、精简 | `Beyond Research`；删除商业书 PDF；精选原创视频 |
| `archive-layout-with-content.md` | 主题演示 | 删除/不发布 | 无 |
| `buildings_halls_in_cuhk.md` | 中文楼宇缩写 | 合并或归档 | `/zh/guides/cuhk-campus/`，写复核日期 |
| `category-archive.html` | 博客分类页 | 按需保留 | 只有 Writing 内容达到一定数量才保留 |
| `collection-archive.html` | 自动列出 collection | 不发布 | 避免暴露未治理 collection |
| `cv.md` | 网页 CV | 全面重写 | 精简结构、唯一数据源、无敏感证明扫描件 |
| `go_to_cuhk.md` | 中文校园出行指南 | 合并或归档 | 与楼宇/到校信息合成一个官方链接优先的指南 |
| `hong_kong_before_arrival.md` | 新生指南 | 归档或重写 | 与 hkinfo 去重；高时效字段必须定期复核 |
| `markdown.md` | Markdown 主题演示 | 删除/不发布 | 无 |
| `mypublication.html` | 主论文页 | 停止手写 | 数据迁入唯一 publication source；旧 URL redirect |
| `non-menu-page.md` | 模板演示 | 删除/不发布 | 无 |
| `page-archive.html` | 自动页面档案 | 不发布 | 使用人工策划 sitemap |
| `portfolio.html` | Portfolio 入口 | 暂停发布 | 权利审查后重建为 Resources |
| `posts.html` | Blog Posts 列表 | 改名 | Writing，按 Research/Personal/Guides 分区 |
| `publications.md` | 第二论文页 | 作为唯一 canonical 页重建 | Selected + filters + complete list |
| `sitemap.md` | 遍历全部内容 | 重写 | 只列经过策划的公开页面 |
| `tag-archive.html` | 标签档案 | 按需保留 | 标签需受控，避免空/低价值页面 |
| `talkmap.html` | 公开假 talk map | 删除/不发布 | 如未来需要，只使用真实 talks 数据重建 |
| `talks.md` | 活动列表 | 重写 | Invited / Oral / Poster / Panels / Attendance |
| `teaching.md` | 课程列表 | 重写 | 角色、责任、主题、成果、合法材料 |
| `terms.md` | 2016 模板隐私条款 | 完全重写 | 根据实际 hosting、analytics、maps、embeds 编写 |
| `travel.md` | Shandong 旅行页 | 与博文合并 | 只留一个 canonical photo essay |

### 30.2 `_posts`

| 文件 | 建议处置 | 具体修改 |
|---|---|---|
| `2023-08-04-blog-post-1.md` | 与 `travel.md` 二选一 | 修正 `/files/trip/` 大小写、补背景/alt/caption、去 EXIF |
| `2023-08-08-hkinfo.md` | 归档或与三个中文 guide 合并 | 删除过期流程和产品背书；只留官方来源；标 `Last verified` |
| `2023-11-16-cuhk_apply.md` | 每年维护或归档 | 更新费用与步骤；修复 HKPFS 链接；删除“随便填”；不推荐第三方处理敏感材料 |
| `2025-07-30-blog-post-1.md` | 保留但深度编辑 | 去重复/刻板印象/Caltech；祖源写 family oral history；加入具体行动；修 permalink |

### 30.3 Collections

| Collection | 当前状态 | 建议 |
|---|---|---|
| `_publications` | 13 条、1 条 YAML 解析异常、缺大量新成果 | 全量重建为唯一数据源，逐 DOI 核验 |
| `_teaching` | 16 条，少 4 门近期课程；日期字段人为占位 | 统一 schema，补 role/responsibilities，删除不安全 syllabus |
| `_portfolio` | 三个高风险文件目录入口 | 立即 `output: false` 或移出发布；权利审查后重建原创 Resources |
| `_talks` | 四条 AcademicPages 假数据 | 删除；未来用真实结构化 talks 数据重新启用 |

### 30.4 数据、模板和公共组件

| 位置 | 建议 |
|---|---|
| `_config.yml` | 更新 title/description/timezone/future/social/analytics；缩小默认 author profile 和公开 collections |
| `_data/navigation.yml` | 改为 7 项核心导航 |
| `_data/authors.yml` | 删除虚构作者；如无多作者需求则移除 |
| `_includes/author-profile.html` | 删除无功能 Follow、ClustrMaps、废弃 HTML 和过量社交链接 |
| `_includes/footer.html` | 模块化并与所有页面统一；只保留核心入口和真实联系方式 |
| `_includes/head/custom.html` | 修复 favicon/manifest，按需加载 MathJax，更新 theme color |
| `_includes/seo.html` | 输出标准 description、Open Graph、Twitter Card、Person JSON-LD |
| `talkmap/` | 删除 debug/fake data，或以后从真实 talk source 生成 |
| `markdown_generator/` | 作为开发工具从生产构建和 sitemap 排除 |
| `_data/comments/` | 删除模板评论；未启用 comments 时不发布相关数据 |

---

## 31. 已知链接与资源修复清单

### 31.1 站内明确失效

- `_posts/2023-08-04-blog-post-1.md` 的 9 张 `/files/trip/...` 图片：实际目录为 `/files/Trip/...`。
- `CUHK_Map_Detail.pdf`：被引用两次，但目标文件缺失；优先改为 CUHK 官方地图链接。
- `_pages/talks.md` 中 `20230708_ustc.pdf`：目标缺失；没有授权/备份就删除链接。
- `_portfolio/Research_Literature_Review.md` 中 `zolfagharian2020closed.pdf`：目标缺失；改 DOI，不补传未知来源 PDF。
- `_publications/zhang2024high2.md`：PDF 指向 `zhai2022survey.pdf`，必须修正。

### 31.2 外部明确失效或低价值

- Undergraduate Course Materials 的旧域名链接中至少 18 个 404；该页面应整体下线，而不是逐个修补第三方课程资料。
- Talk map 上游 notebook 链接 404；页面本身应删除。
- Talks 中一个研讨会页面 404；用官方存档、DOI 或删除链接。
- hkinfo 中旧校巴应用 404；用 CUHK 官方交通页并写复核日期。

### 31.3 图标和 manifest

- `images/manifest.json` 声明的 6 个 Android 图标不存在。
- `_includes/head/custom.html` 引用的 `favicon.ico` 不存在。
- manifest 名称仍为 `Minimal Mistakes`，应改为 `Liuchao Jin | Homepage`。

### 31.4 链接迁移规则

1. 站内文件使用 Jekyll `relative_url` 或根相对路径。
2. 文件夹和文件名统一小写短横线，避免 Windows 可用、Linux 404。
3. 当前站点的 507 个 HTTP 绝对链接逐步改为站内相对路径。
4. 其他 81 个 HTTP 外链优先升级 HTTPS；没有 HTTPS 且不再维护的链接替换或删除。
5. 旧域名统一 redirect，不在正文和 PDF 中继续硬编码。
