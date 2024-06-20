import numpy as np

# Добавил структуру, пытаясь придумать как хранить исходные данные и  
# промежуточные профили выравнивания без лишнего пересчета
class Sequence:
    """ Примитивный pfm на int8 для хранения данных """
    def __init__(self, seq, number):
        # Вызов базового конструктора
        if (seq == ""):
            self.len = number
            self.data = np.zeros((5, number), dtype=np.int8)
            self.is_consensus = False
            self.seq_count = 0
            self.parents = []

        # если в конструктор передана какая-то последовательность
        else:
            self.len = len(seq)
            self.data = np.zeros((5, self.len), dtype=np.int8)
            for i in range(self.len):
                if seq[i] == "A":
                    self.data[0][i] += 1
                elif seq[i] == "C":
                    self.data[1][i] += 1
                elif seq[i] == "G":
                    self.data[2][i] += 1
                elif seq[i] == "T":
                    self.data[3][i] += 1
                elif seq[i] == "-":
                    self.data[4][i] += 1

            self.is_consensus = False
            self.seq_count = 1
            self.parents = [number]

    def get_max_freq(self, col):
        """
        Вспомогательная функция для подсчета редакционного расстояния.
        Возвращает индексы тех строк, которые содержат наибольшую частоту в заданном столбце матрицы
        """
        max_frequency = np.max(self.data[:][col])
        return np.where(self.data[:][col] == max_frequency)

    def concat_freq(self, freq, col):
        """ Добавляет столбец freq к значениям столбца col """
        self.data[:, col] += freq

    def add_gap(self, col):
        """ Добавляет гэп в позицию col """
        self.data[4][col] += 1

    def update_fields(self, seq1, seq2):
        """ Внесем информацию для консенсусной последовательности, если вызвали базовый конструктор"""
        self.is_consensus = True
        self.seq_count = seq1.seq_count + seq2.seq_count
        self.parents = seq1.parents + seq2.parents



    # def merge_with_seq(self, seq):
    #     self.is_consensus = True
    #     self.seq_count += seq.seq_count
    #     self.parents.append(seq.seq_count)

    #     if (self.data.shape[1] == seq.data.shape[1]):
    #         self.data += seq.data

    # # с осторожностью, если объединяем консенсусы
    # def merge_with_consensus(self, seq):
    #     self.is_consensus = True
    #     self.seq_count += seq.seq_count
    #     self.parents.append(seq.seq_count)

    #     if (self.data.shape[1] == seq.data.shape[1]):
    #         self.data += seq.data
    #     # если можно дописать в текущую матрицу
    #     elif (self.data.shape[1] > seq.data.shape[1]):
    #         for i in range(seq.data.shape[1]):
    #             self.data[:][i] += seq.data[:][i]
    #     # нужно увеличивать число столбцов матрицы
    #     else:
    #         columns = np.zeros((4, seq.data.shape[1] - self.data.shape[1])) 
    #         self.data = np.hstack((self.data, columns))
    #         for i in range(seq.data.shape[1]):
    #             self.data[:][i] += seq.data[:][i]

    #     ####
    #     ####
    #     # нужно подумать что делать если объединили 2 консенсуса с общим родителем





class MultAlign:
    def __init__(self, seq_list,
                 delete_penalty = 1,
                 insert_penalty = 1,
                 mismatch_penalty = 1,
                 match_reward = 0):

        self.delete_penalty = delete_penalty,
        self.insert_penalty = insert_penalty,
        self.mismatch_penalty = mismatch_penalty,
        self.match_reward = match_reward 

        self.sequences = [Sequence(seq_list[i], i) for i in range(len(seq_list))]


    # К 1 дз было замечание: скрипт должен потреблять O(min(n, m)) памяти
    # Учел
    def vagner_fischer(self, seq1, seq2):
        """ Функция для подсчета редакционного расстояния. Принимает строки, штрафы и награду за совпадение """
        n, m = seq1.len, seq2.len

        # При подсчете будем считать, что первая последовательность длинее
        # т.е. матрица dp вытянута вниз
        if (n < m):
            n, m = m, n
            seq1, seq2 = seq2, seq1

        # Инициализация
        dp = np.zeros((2, m + 1), dtype=int)
        dp[1][0] = 1

        for j in range(m + 1):
            dp[0][j] = j

        for i in range(1, n + 1):
            # на i-й итерации цикла мы рассматриваем i-1 и i строки исходной матрицы n*m  
            for j in range(1, m + 1):
                # если был подан консенсус с коллизией, то проверяем наличие общего нуклеотида - далее коллизия будет устранена при слиянии
                match_count = np.intersect1d(seq1.get_max_freq(i-1), seq2.get_max_freq(j-1)).size

                score = self.match_reward if ( match_count > 0 ) else self.mismatch_penalty
                dp[1][j] = min(dp[0][j] + self.delete_penalty,
                            dp[1][j - 1] + self.insert_penalty,
                            dp[0][j - 1] + score)

            # стираем более ненужную i-1 строку исходной матрицы
            # и задаем инициализирующее значение для i-й строки
            dp = dp[::-1, ::]
            dp[1][0] = i

        return dp[1][m]


    # так и не понял какое полагается поведение при выравнивании консенсуса с коллизией
    # поэтому не стал пытаться усложнять и реализовал классический алгоритм
    def needleman_wunsch(self, seq1, seq2,
                         gap_penalty = -2,
                         mismatch_penalty = -1,
                         match_reward = 1
                         ):
        """
        Алгоритм Нидлмана-Вунша для получения консенсусной последовательности,
        хранящейся в матрице частот
        """
        align = np.zeros((seq1.len + 1, seq2.len + 1), dtype=int)
        align[:, 0] = np.arange(0, gap_penalty*align.shape[0], gap_penalty)
        align[0, :] = np.arange(0, gap_penalty*align.shape[1], gap_penalty)
        
        # Для отслеживания будем записывать ось по которой нужно двигаться, восстанавливая путь
        # 0 - движение по вертикали  <==>  гэп у seq2
        # 1 - по горизонтали  <==>  гэп у seq1
        # 2 - по диагонали
        traceback_matrix = np.zeros((seq1.len + 1, seq2.len + 1), dtype=int)
        traceback_matrix[0][:] = np.ones((1, traceback_matrix.shape[1]))
        
        for i in range(1, align.shape[0]):
            for j in range(1, align.shape[1]):

                if (np.intersect1d(seq1.get_max_freq(i-1), seq2.get_max_freq(j-1)).size > 0):
                    diag_score = align[i-1][j-1] + match_reward
                else:
                    diag_score = align[i-1][j-1] + mismatch_penalty

                up_score = align[i-1][j] + gap_penalty
                left_score = align[i][j-1] + gap_penalty

                best_score = max(diag_score, up_score, left_score)
                if (best_score == diag_score):
                    align[i][j] = best_score
                    traceback_matrix[i][j] = 2

                elif (best_score == left_score):
                    align[i][j] = best_score
                    traceback_matrix[i][j] = 1

                elif (best_score == up_score):
                    align[i][j] = best_score
                    traceback_matrix[i][j] = 0



        # восстановление пути
        # q1: как инициализировать матрицу консенсусной строки
        # q2: как вообще считать в случае с матрицами

        #  сначала уточним обратный маршрут
        traceback_way = []
        while (i > 0) or (j > 0):
            if (i > 0) and (j > 0) and (traceback_matrix[i][j] == 2):
                traceback_way.append(2)
                i -= 1
                j -= 1
            elif (i > 0) and (traceback_matrix[i][j] == 0):
                traceback_way.append(0)
                i -= 1
            else:
                traceback_way.append(1)
                j -= 1
        del traceback_matrix
        traceback_way.reverse()

        # затем попробуем собрать матрицу, включая в неё столбцы из матриц исходных последовательностей
        # при необходимости будем инкрементировать ячейки с гэпами и пропускать итерацию цикла на соответствующей последовательности
        consensus = Sequence(seq="", number=len(traceback_way))

        # при восстановлении последовательности
        # 0 - движение по вертикали  <==> добавляем гэп к seq2
        # 1 - движение по горизонтали  <==>  добавляем гэп к seq1

        seq1_gaps_count, seq2_gaps_count = 0, 0
        for i in range(1, len(traceback_way)+1):
            way = traceback_way.pop()
            if (way == 2):
                consensus.concat_freq(seq1.len - i + seq1_gaps_count)
                consensus.concat_freq(seq2.len - i + seq2_gaps_count)
            elif (way == 1):
                # сдвиг по горизонтали - у первой последовательности гэп
                consensus.add_gap(seq1.len - i + seq1_gaps_count)
                seq1_gaps_count += 1

                consensus.concat_freq(seq2.len - i + seq2_gaps_count)                
            else:
                # сдвиг по вертикали - у второй последовательности гэп
                consensus.concat_freq(seq1.len - i + seq1_gaps_count)                

                consensus.add_gap(seq2.len - i + seq2_gaps_count)
                seq2_gaps_count += 1

        consensus.update_fields(seq1, seq2)
        return consensus


    def find_closest_pair(self, strings_array):
        """ Функция для нахождения двух самых близких строк в массиве """
        min_distance, closest_pair = float('inf'), None

        # Проходим по всем парам строк и ищем минимальное расстояние левенштейна
        for i in range(len(strings_array)):
            for j in range(i + 1, len(strings_array)):
                distance = self.vagner_fischer(strings_array[i], strings_array[j])
                if distance < min_distance:
                    min_distance = distance
                    closest_pair = (i, j)
        return closest_pair


    def generate_consensus_string(self, strings):
        """ Функция для создания консенсусной строки из массива строк """
        consensus = ''

        # Проходим по символам в каждой позиции и выбираем наиболее часто встречающийся символ
        for i in range(len(strings[-1])):
            counts = {}
            for string in strings:
                if string[i] not in counts:
                    counts[string[i]] = 1
                else:
                    counts[string[i]] += 1
            consensus += max(counts, key=counts.get)

        return consensus


    def greedy_mult_align(self):
        """ Жадное множественное выравнивание """
        alignments = []
        strings = self.sequences
        
        while len(strings) > 1:
            # Находим ближайшую пару строк
            i, j = self.find_closest_pair(strings)

            # Гененируем консенсусную строку для пары i,j
            # consensus = self.generate_consensus_string([strings[i], strings[j]])
            consensus = self.needleman_wunsch(strings[i], strings[j])

            # Добавляем выравнивание в список
            ### alignments.append((strings[i], strings[j], consensus))
            alignments.append(consensus)
            
            # Удаляем использованные строки из массива
            del strings[j]
            del strings[i]

            # Добавляем консенсусную строку в массив
            strings.append(consensus)

        # Возвращаем результат выравнивания
        return alignments




###############################################################################
# Пример из лекции
strings = ["GATTCA", "GTCTGA", "GATATT", "GTCAGC"]

delete_penalty = 1
insert_penalty = 1
mismatch_penalty = 1
match_reward = -1

needleman_wunsch_match_reward = 1
needleman_wunsch_mismatch_penalty = -1
needleman_wunsch_gap_penalty = -2

alignments = MultAlign(strings)
alignments.greedy_mult_align()
needleman_wunsch_gap_penalty()
